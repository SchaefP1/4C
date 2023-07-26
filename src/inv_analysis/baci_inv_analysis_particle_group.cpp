/*----------------------------------------------------------------------*/
/*! \file
\brief Collection of particles for the Sequential Monte Carlo algorithm

\level 3

*/
/*----------------------------------------------------------------------*/

#include "baci_inv_analysis_particle_group.H"

#include "baci_inv_analysis_base.H"
#include "baci_inv_analysis_utils.H"
#include "baci_inv_analysis_likelihood_evaluation.H"
#include "baci_inv_analysis_metropolis_kernel.H"
#include "baci_inv_analysis_particle_comm.H"
#include "baci_inv_analysis_particle_data.H"
#include "baci_lib_globalproblem.H"
#include "baci_linalg_utils_sparse_algebra_math.H"

#include <iomanip>
#include <list>
#include <random>
#include <stdexcept>

#ifdef TRAP_FE
#include <cfenv>
#endif


/*----------------------------------------------------------------------*/
INVANA::ParticleGroup::ParticleGroup(const Teuchos::ParameterList& invp)
    : particle_data_(),
      weights_(),
      loglikemixture_(Teuchos::null),
      lhcorrectpost_(0.0),
      lhcorrectprior_(0.0),
      mc_kernel_(Teuchos::null),
      mc_adapt_scale_(0.1),
      mc_kernel_iter_(1),
      gnumparticles_(0),
      lnumparticles_(0),
      ngroups_(0),
      mygroup_(0),
      params_(invp)
{
}

/*----------------------------------------------------------------------*/
void INVANA::ParticleGroup::Init(Teuchos::RCP<LogLikeMixture> mixtures)
{
  loglikemixture_ = mixtures;

  // init mc kernel
  mc_kernel_ = Teuchos::rcp(new MetropolisKernel());
  mc_kernel_->Init(mixtures);
  mc_adapt_scale_ = params_.get<double>("MC_INIT_SCALE");
  mc_kernel_iter_ = params_.get<int>("SMC_KERNEL_ITER");

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::ParticleGroup::Setup()
{
  DRT::Problem* problem = DRT::Problem::Instance();

  // get some global numbers and validate
  ngroups_ = problem->GetCommunicators()->NumGroups();
  mygroup_ = problem->GetCommunicators()->GroupId();
  gnumparticles_ = params_.get<int>("NUM_PARTICLES");

  // only allow for euqally distributed particles among groups
  if (gnumparticles_ % ngroups_)
    dserror("Choose ngroups % nparticles == 0 to allow for euqual distribution");

  // local number of particles in each group
  lnumparticles_ = (int)gnumparticles_ / ngroups_;

  // setup communicators
  SetupComms();

  // assign gids to particles
  InitializeParticleData();

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::ParticleGroup::SetupComms()
{
  DRT::Problem* problem = DRT::Problem::Instance();

  Teuchos::RCP<Epetra_Comm> gcomm = problem->GetCommunicators()->GlobalComm();
  Teuchos::RCP<Epetra_Comm> lcomm = problem->GetCommunicators()->LocalComm();

  // Construct inter group communicator
  pcomm_ = Teuchos::rcp(new ParticleComm());
  pcomm_->Init(gcomm, lcomm, lnumparticles_, mygroup_);
  pcomm_->Setup();

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::ParticleGroup::InitializeParticleData()
{
  // particles within each group are offset by the number of groups
  for (int i = 0; i < lnumparticles_; i++) my_particle_gids_.push_back(mygroup_ + i * ngroups_);

  // set up empty data
  for (int i = 0; i < lnumparticles_; i++)
  {
    int pgid = my_particle_gids_[i];
    weights_[pgid] = 1.0 / gnumparticles_;
    Data()[pgid] = Teuchos::rcp(new INVANA::ParticleData());
    Data()[pgid]->Init(loglikemixture_->StateMap());
  }
  new_weights_ = weights_;

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::ParticleGroup::PreEvaluate()
{
  // no need for particular fpe trapping in here; if the particle
  // evaluation at the MAP-point fails something else is wrong!

  // Get mean from prior as preevaluation point
  Epetra_Vector sample(loglikemixture_->EvalPrior().Mean());

  double posterior;
  double prior;
  int err = 0;

  // evaluate mixture loglikelihood at the sample
  err = EvaluateMixture(sample, posterior, prior);
  if (err) dserror("PreEvaluation of particles failed. This is fatal!");

  // compute likelihood correction factor
  lhcorrectpost_ = posterior;
  lhcorrectprior_ = prior;

  // lets wait here
  pcomm_->GComm().Barrier();

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::ParticleGroup::DrawInitialStates()
{
#ifdef TRAP_FE
  fedisableexcept(FE_ALL_EXCEPT);
#endif

  // Draw intial states from prior
  Epetra_Vector sample(loglikemixture_->StateMap(), false);
  for (DATAITER it = Data().begin(); it != Data().end(); it++)
  {
    double posterior;
    double prior;
    int iraised = 0;
    int graised = 1;
    int err = 0;
    while (graised)
    {
      try
      {
        // clear exception being already signaled
#ifdef TRAP_FE
        feclearexcept(FE_ALL_EXCEPT);
#endif

        // get new proposal state
        loglikemixture_->DrawfromPrior(sample);

        // evaluate mixture loglikelihood at the sample
        err = EvaluateMixture(sample, posterior, prior);

#ifdef TRAP_FE
        // test for occurence of these signals
        iraised = fetestexcept(FE_INVALID | FE_DIVBYZERO);
#endif

        // and make them known in the local world
        sample.Comm().SumAll(&iraised, &graised, 1);

        // throw this error now consistently across the local group
        if (graised or err)
          throw std::runtime_error("FPE/Convergence failure during mixture evaluation");
      }
      catch (std::exception& e)
      {
        std::cout << "(Group " << mygroup_
                  << ") Baci was not able to compute this sample due to:" << std::endl;
        std::cout << e.what() << std::endl;
        std::cout << "-> draw another sample" << std::endl;
        graised = 1;
        sample.Comm().Barrier();
      }
    }
    // if successfull set to ParticleData
    it->second->SetState(sample);
    it->second->SetData(posterior, prior);
  }

  // lets wait here
  pcomm_->GComm().Barrier();

#ifdef TRAP_FE
  feclearexcept(FE_ALL_EXCEPT);
  feenableexcept(FE_INVALID | FE_DIVBYZERO);
#endif

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::ParticleGroup::SetData(DATA& data)
{
  for (DATAITER it = Data().begin(); it != Data().end(); it++)
  {
    if (not data.count(it->first)) dserror("key mismatch!");
    *(it->second) = *data[it->first];
  }

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::ParticleGroup::SetWeights(std::map<int, double>& weights)
{
  weights_ = weights;
  new_weights_ = weights;

  return;
}

/*----------------------------------------------------------------------*/
const INVANA::DATA& INVANA::ParticleGroup::GetData() { return Data(); }

/*----------------------------------------------------------------------*/
const std::map<int, double>& INVANA::ParticleGroup::GetWeights() { return weights_; }

/*----------------------------------------------------------------------*/
int INVANA::ParticleGroup::EvaluateMixture(
    const Epetra_Vector& state, double& posterior, double& prior)
{
  int err = loglikemixture_->EvaluateMixture(state, posterior, prior);
  return err;
}

/*----------------------------------------------------------------------*/
double INVANA::ParticleGroup::NewEffectiveSampleSize(double scale_next, double scale_curr)
{
  // in here i dont want SIGFPE but a more sensible treatment
#ifdef TRAP_FE
  fedisableexcept(FE_ALL_EXCEPT);
  feclearexcept(FE_ALL_EXCEPT);
#endif
  int raise = 0;
  int graise = 0;

  // --------------- Compute new weights
  // unnormalized weights
  std::vector<double> myweights;
  // and their ids
  std::vector<int> myids;

  // loop all the particles in this group
  for (DATAITER it = Data().begin(); it != Data().end(); it++)
  {
    // get posterior value
    double posterior = it->second->GetPosterior();
    posterior -= lhcorrectpost_;

    // get prior value
    double prior = it->second->GetPrior();
    prior -= lhcorrectprior_;

    // compute current mixture
    double m_curr = posterior * scale_curr + prior * (1.0 - scale_curr);
    // compute next mixture
    double m_next = posterior * scale_next + prior * (1.0 - scale_next);

    // compute likelihood ratio
    double val = exp(m_next - m_curr);  // possibly overflows
    myweights.push_back(weights_[it->first] * val);

    // keep track of the order of id (just out of paranoia)
    myids.push_back(it->first);
  }

  // normalize
  std::vector<double> weights(gnumparticles_);
  pcomm_->IComm().GatherAll(myweights.data(), weights.data(), lnumparticles_);

  double sum = 0.0;
  for (int i = 0; i < (int)weights.size(); i++) sum += weights[i];

  // set to new_weights_
  int i = 0;
  for (std::vector<int>::iterator it = myids.begin(); it != myids.end(); it++)
  {
    new_weights_[*it] = myweights[i] / sum;
    i++;
  }
  // --------------- END Compute new weights

  //  --------------- Compute ESS
  // compute the same ess on all procs
  double ess = 0.0;
  double sum1 = 0.0;
  double sum2 = 0.0;
  for (int i = 0; i < (int)weights.size(); i++)
  {
    sum1 += weights[i];
    sum2 += weights[i] * weights[i];
  }
  // sum1 is actually one since we have normalized weights
  ess = sum1 * sum1 / sum2;
  //  --------------- END Compute ESS

  // test flags
#ifdef TRAP_FE
  raise = fetestexcept(FE_UNDERFLOW | FE_OVERFLOW | FE_DIVBYZERO);
#endif
  pcomm_->GComm().SumAll(&raise, &graise, 1);

  // suppose the step was just too large! If it
  // wasn't we can not use this computer anyways!
  if (graise)
  {
    if (pcomm_->GComm().MyPID() == 0)
      std::cout << "caught invalid operation setting ess = 0.0" << std::endl;

    ess = 0.0;
  }

  // reset flags and enable traps in case
#ifdef TRAP_FE
  feclearexcept(FE_ALL_EXCEPT);
  feenableexcept(FE_INVALID | FE_DIVBYZERO);
#endif

  // everybody should be done before proceeding
  pcomm_->GComm().Barrier();

  return ess;
}

/*----------------------------------------------------------------------*/
void INVANA::ParticleGroup::UpdateWeights()
{
  weights_ = new_weights_;

  // every proc must have updated before proceeding
  pcomm_->GComm().Barrier();

  return;
}


/*----------------------------------------------------------------------*/
double INVANA::ParticleGroup::EffectiveSampleSize()
{
  // no need for switching off the fpe traps, since only
  // numerically valid sets of weights should end up in
  // weights_ anyways

  double ess = 0.0;

  // put weights into vector
  std::vector<double> myweights;
  std::map<int, double>::iterator it;
  for (it = weights_.begin(); it != weights_.end(); it++) myweights.push_back(it->second);

  // gather all weights on all procs
  std::vector<double> weights(gnumparticles_);
  pcomm_->IComm().GatherAll(myweights.data(), weights.data(), lnumparticles_);

  // compute the same ess on all procs
  double sum = 0.0;
  double sum2 = 0.0;
  for (int i = 0; i < (int)weights.size(); i++)
  {
    sum += weights[i];
    sum2 += weights[i] * weights[i];
  }
  // sum is actually one since we have normalized weights
  ess = sum * sum / sum2;

  return ess;
}

/*----------------------------------------------------------------------*/
void INVANA::ParticleGroup::ResampleParticles()
{
  // get all particle ids and weights
  std::vector<double> myweights;
  for (int i = 0; i < lnumparticles_; i++) myweights.push_back(weights_[my_particle_gids_[i]]);

  // all particles and weights an all procs
  std::vector<int> particles(gnumparticles_);
  std::vector<double> weights(gnumparticles_);
  pcomm_->IComm().GatherAll(my_particle_gids_.data(), particles.data(), lnumparticles_);
  pcomm_->IComm().GatherAll(myweights.data(), weights.data(), lnumparticles_);

  std::vector<double> intervals(gnumparticles_ + 1);
  for (int i = 0; i <= gnumparticles_; i++) intervals[i] = (double)i;

  // generate the distribution
  std::default_random_engine generator;
  std::piecewise_constant_distribution<double> distribution(
      intervals.begin(), intervals.end(), weights.begin());

  // draw gnumparticles
  std::vector<int> draw(gnumparticles_);
  for (int i = 0; i < gnumparticles_; i++) draw[i] = (int)distribution(generator);

  // get the particle gids which survived
  for (int i = 0; i < gnumparticles_; i++) particles[i] = draw[i];

  RedistributeParticleData(particles);

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::ParticleGroup::RedistributeParticleData(std::vector<int> pgids)
{
  //  if (pcomm_->GComm().MyPID()==0)
  //  {
  //    std::cout << "Particles to be redistributed: ";
  //    for (int i=0; i<(int)pgids.size(); i++)
  //      std::cout << pgids[i] << " ";
  //    std::cout << std::endl;
  //  }

  if ((int)pgids.size() != gnumparticles_)
    dserror("something went wrong in resampling the particles!");

  // reorganize pgids as a map <pgids, # to be sent>
  std::map<int, int> pcount;
  for (int i = 0; i < gnumparticles_; i++)
  {
    if (pcount.count(pgids[i]))
      pcount[pgids[i]] += 1;
    else
      // count starts with 0 since the sender already as one
      pcount[pgids[i]] = 0;
  }

  std::map<int, int>::iterator it;

  // promote free spots to be exported to
  // initialize with the potentially lnumparticles_ free spots
  std::list<int> nrecvs(my_particle_gids_.begin(), my_particle_gids_.end());
  for (it = pcount.begin(); it != pcount.end(); it++) nrecvs.remove(it->first);

  //  std::cout << "mylist contains:";
  //  for (auto it=nrecvs.begin(); it!=nrecvs.end(); ++it)
  //    std::cout << ' ' << *it;
  //  std::cout << '\n';

  // loop the particles to be distributed (every group and proc knows them)
  DATA particle_data;
  for (it = pcount.begin(); it != pcount.end(); it++)
  {
    int ptosend = it->first;
    int ntosend = it->second;

    // if there is only one -> just keep it
    if (ntosend == 0) continue;

    // every group having a free spot should take one or more
    int irecv = 0;
    std::vector<int> pgidtarget;
    int ntorecv = 0;
    int stilltosend = ntosend;
    for (int j = 0; j < ngroups_; j++)
    {
      if (mygroup_ == j)
      {
        // if there are free spots in this group
        if (nrecvs.size() > 0)
        {
          int nn = 0;
          if (stilltosend <= (int)nrecvs.size())
            nn = stilltosend;  // this group takes all
          else if (stilltosend > (int)nrecvs.size())
            nn = (int)nrecvs.size();  // this group takes some
          else
            dserror("don't know this case!");

          irecv = 1;
          for (int i = 0; i < nn; i++)
          {
            pgidtarget.push_back(nrecvs.back());
            nrecvs.pop_back();
          }

          // increase number of received particles
          ntorecv += nn;
          stilltosend -= nn;

        }  // nrecvs>0
      }    // mygroup

      // broadcast number of received particles
      pcomm_->IComm().Broadcast(&ntorecv, 1, j);
      pcomm_->IComm().Broadcast(&stilltosend, 1, j);
      // pcomm_->GComm().Barrier();

      // all particles have a destination already skip the rest of the groups
      if (ntorecv == ntosend) break;
    }

    // build from map and set data
    int isend = 0;
    std::map<int, Teuchos::RCP<ParticleData>> data;
    int dummygid = 0;
    double post = 0.0;
    double prior = 0.0;
    double state = 0.0;
    if (weights_.count(ptosend))
    {
      isend = 1;
      // we can just set the pointer here!
      // copying is done upon communication
      data[dummygid] = Data()[ptosend];

      // set data to compare
      post = data[0]->GetPosterior();
      prior = data[0]->GetPrior();
      data[0]->GetState().Norm2(&state);
    }
    Epetra_Map frommap(-1, isend, &dummygid, 0, pcomm_->IComm());

    // build the tomap
    Epetra_Map tomap(-1, irecv, &dummygid, 0, pcomm_->IComm());

    // export
    DRT::Exporter ex(frommap, tomap, pcomm_->IComm());
    ex.Export(data);

    // broadcast particle data to compare against
    double postall = 0.0;
    double priorall = 0.0;
    double stateall = 0.0;
    pcomm_->IComm().SumAll(&post, &postall, 1);
    pcomm_->IComm().SumAll(&prior, &priorall, 1);
    pcomm_->IComm().SumAll(&state, &stateall, 1);

    // store at the new positions if i have any new particles
    for (std::vector<int>::iterator kt = pgidtarget.begin(); kt != pgidtarget.end(); kt++)
    {
      particle_data[*kt] = data[dummygid];

      double priornew = particle_data[*kt]->GetPrior();
      double postnew = particle_data[*kt]->GetPosterior();
      double normstate;
      particle_data[*kt]->GetState().Norm2(&normstate);
      double eps = 1.0e-14;
      if (abs(priorall - priornew) > eps)
      {
        std::cout << std::setprecision(10) << "priorall: " << priorall << " priornew: " << priornew
                  << std::endl;
        dserror("Data distribution failed: Wrong prior data.");
      }
      if (abs(postall - postnew) > eps) dserror("Data distribution failed: Wrong posterior data.");
      if (abs(stateall - normstate) > eps) dserror("Data distribution failed: Wrong state data.");
    }
    pcomm_->GComm().Barrier();
  }

  // update ParticleData where necessary
  for (DATAITER it = particle_data.begin(); it != particle_data.end(); it++)
    Data()[it->first] = it->second;

  // reinitialize weights
  for (std::map<int, double>::iterator it = weights_.begin(); it != weights_.end(); it++)
    it->second = 1.0 / gnumparticles_;
}

/*----------------------------------------------------------------------*/
void INVANA::ParticleGroup::RejuvenateParticles(double scale)
{
  int numiter = mc_kernel_iter_;
  double covscale = mc_adapt_scale_;

  // backup copy of the particle data (in case the MC move is not proper)
  std::map<int, Teuchos::RCP<ParticleData>> p_bak;
  for (DATAITER it = Data().begin(); it != Data().end(); it++)
  {
    p_bak[it->first] = Teuchos::rcp(new ParticleData(*(it->second)));
  }

  double boundl = 0.20;
  double boundu = 0.40;
  double acceptance = 0.0;
  while (acceptance < boundl || acceptance > boundu)
  {
    acceptance = 0.0;

    // loop all the particles in this group
    double acc;
    for (DATAITER it = Data().begin(); it != Data().end(); it++)
    {
      mc_kernel_->Sample(numiter, scale, covscale, *(it->second), acc);
      acceptance += acc;
    }

    // compute global acceptance rate;
    double acc_sum;
    pcomm_->IComm().SumAll(&acceptance, &acc_sum, 1);
    acceptance = acc_sum / gnumparticles_;

    if (pcomm_->GComm().MyPID() == 0)
    {
      printf("  MCMC Move: propscal %.3f, acc: %.2f\n", covscale, acceptance);
      fflush(stdout);
    }

    if (acceptance < boundl)
    {
      covscale *= 0.5;
      // reset
      SetData(p_bak);
    }
    else if (acceptance > boundu)
    {
      covscale *= 1.2;
      // reset
      SetData(p_bak);
    }
  }

  mc_adapt_scale_ = covscale;

  return;
}
/*----------------------------------------------------------------------*/
void INVANA::ParticleGroup::ComputeMean(Epetra_Vector& mean, Epetra_Vector& stdev)
{
  std::map<int, Teuchos::RCP<ParticleData>> data = Data();

  ComputeMean(data, mean, stdev);
}
/*----------------------------------------------------------------------*/
void INVANA::ParticleGroup::ComputeMean(
    std::map<int, Teuchos::RCP<ParticleData>>& data, Epetra_Vector& mean, Epetra_Vector& stdev)
{
  // data to be communicated
  std::map<int, double> weights = weights_;

  // particle ids
  std::vector<int> pgids(gnumparticles_);
  pcomm_->IComm().GatherAll(my_particle_gids_.data(), pgids.data(), lnumparticles_);

  // set up tomap and frommap
  int ntosend = lnumparticles_;
  int ntorecv = 0;
  if (mygroup_ == 0) ntorecv = gnumparticles_;

  Epetra_Map frommap(-1, ntosend, my_particle_gids_.data(), 0, pcomm_->IComm());
  Epetra_Map tomap(-1, ntorecv, pgids.data(), 0, pcomm_->IComm());

  // export
  DRT::Exporter ex(frommap, tomap, pcomm_->IComm());
  ex.Export(data);
  ex.Export(weights);

  if (mygroup_ == 0)
  {
    // check for normalization
    double sum = 0.0;
    for (std::map<int, double>::iterator it = weights.begin(); it != weights.end(); it++)
      sum += it->second;

    if (abs(sum - 1.0) > 1.0e-1) dserror("weight is not normalized, but %.1f", sum);

    // mean
    mean.Scale(0.0);
    for (DATAITER it = data.begin(); it != data.end(); it++)
    {
      double w = weights[it->first];
      mean.Update(w, it->second->GetState(), 1.0);
    }

    // standard deviation
    stdev.Scale(0.0);
    Epetra_Vector dummy(stdev.Map(), false);
    for (DATAITER it = data.begin(); it != data.end(); it++)
    {
      // particle-mean
      dummy.Scale(1.0, it->second->GetState());
      dummy.Update(-1.0, mean, 1.0);

      // (particle-mean).^2
      double* val;
      dummy.ExtractView(&val);
      for (int i = 0; i < dummy.MyLength(); i++) val[i] = val[i] * val[i];

      double w = weights[it->first];
      stdev.Update(w, dummy, 1.0);
    }
    // sqrt(sum_i(w_i*(particle_i-mean)^2))
    double* val;
    stdev.ExtractView(&val);
    for (int i = 0; i < dummy.MyLength(); i++) val[i] = sqrt(val[i]);
  }

  return;
}