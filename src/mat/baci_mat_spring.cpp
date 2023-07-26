/*----------------------------------------------------------------------*/
/*! \file
\brief Material law for elastic spring (either translational or rotational spring)

\level 3

*----------------------------------------------------------------------*/


#include <vector>
#include "baci_mat_spring.H"
#include "baci_lib_globalproblem.H"
#include "baci_mat_par_bundle.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::Spring::Spring(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      stiffness_(matdata->GetDouble("STIFFNESS")),
      density_(matdata->GetDouble("DENS"))
{
}


Teuchos::RCP<MAT::Material> MAT::PAR::Spring::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Spring(this));
}

MAT::SpringType MAT::SpringType::instance_;


DRT::ParObject* MAT::SpringType::Create(const std::vector<char>& data)
{
  MAT::Spring* spring = new MAT::Spring();
  spring->Unpack(data);
  return spring;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Spring::Spring() : params_(NULL) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Spring::Spring(MAT::PAR::Spring* params) : params_(params) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Spring::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Spring::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::Spring*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}