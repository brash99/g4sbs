#include "G4SBSICHit.hh"
//______________________________________________________________________________
G4ThreadLocal G4Allocator<G4SBSICHit>* G4SBSICHitAllocator = 0;
//______________________________________________________________________________
G4SBSICHit::G4SBSICHit()
 : G4VHit(),
   fEdep(0.),
   fTrackLength(0.),
   fEtot(0.),
   fBeta(0.),
   fHitTime(0.),
   fTrackID(-1),
   fPID(-1),
   fMID(-1)
{
   fPos.setX(0);    fPos.setY(0);    fPos.setZ(0);
   fLabPos.setX(0); fLabPos.setY(0); fLabPos.setZ(0);
   fMom.setX(0);    fMom.setY(0);    fMom.setZ(0);
   fPrintToCSV = false;
   fCntr = 0;
}
//______________________________________________________________________________
G4SBSICHit::~G4SBSICHit()
{

}
//______________________________________________________________________________
G4SBSICHit::G4SBSICHit(const G4SBSICHit& right)
  : G4VHit()
{
  fEdep        = right.fEdep;
  fTrackLength = right.fTrackLength;
  fEtot        = right.fEtot;
  fBeta        = right.fBeta;
  fHitTime     = right.fHitTime;
  fPID         = right.fPID;
  fMID         = right.fMID;
  fPos         = right.fPos;
  fLabPos      = right.fLabPos;
  fMom         = right.fMom;
  fPrintToCSV  = right.fPrintToCSV;
  fCntr        = right.fCntr; 
}
//______________________________________________________________________________
const G4SBSICHit& G4SBSICHit::operator=(const G4SBSICHit& right)
{
  fEdep        = right.fEdep;
  fTrackLength = right.fTrackLength;
  fEtot        = right.fEtot;
  fBeta        = right.fBeta;
  fHitTime     = right.fHitTime;
  fPID         = right.fPID;
  fMID         = right.fMID;
  fPos         = right.fPos;
  fLabPos      = right.fLabPos;
  fMom         = right.fMom;
  fPrintToCSV  = right.fPrintToCSV;
  fCntr        = right.fCntr; 
  return *this;
}
//______________________________________________________________________________
G4bool G4SBSICHit::operator==(const G4SBSICHit& right) const
{
  return ( this == &right ) ? true : false;
}
//______________________________________________________________________________
void G4SBSICHit::Print()
{
  G4cout << "Edep: "
         << std::setw(7) << G4BestUnit(fEdep,"Energy")
         << " track length: "
         << std::setw(7) << G4BestUnit(fTrackLength,"Length")
         << " track x: "
         << std::setw(7) << G4BestUnit(fPos.getX(),"Length")
         << " track y: "
         << std::setw(7) << G4BestUnit(fPos.getY(),"Length")
         << " track z: "
         << std::setw(7) << G4BestUnit(fPos.getZ(),"Length") << G4endl;

  if(fPrintToCSV) PrintToCSV(); 
}
//______________________________________________________________________________
void G4SBSICHit::PrintToCSV(){

  char outpath[200],msg[200]; 
  sprintf(outpath,"ionChamber_hits.csv");

  std::ofstream outfile;
  outfile.open(outpath); 
  if(outfile.fail() ){
     if(fCntr==1) std::cout << "[G4SBSICHit]: Cannot open the file " << outpath << std::endl;
     fCntr++; 
  }else{
     // if(fCntr==0){
     //    sprintf(msg,"Edep_keV,trackLen_mm,trackX_mm,trackY_mm,trackZ_mm"); 
     //    outfile << msg << std::endl; 
     // } 
     sprintf(msg,"%.3lf,%.3lf,%.3lf,%.3lf,%.3lf",
	     fEdep/CLHEP::keV,fTrackLength/CLHEP::mm,fPos.getX()/CLHEP::mm,fPos.getY()/CLHEP::mm,fPos.getZ()/CLHEP::mm);
     outfile << msg << std::endl;
     outfile.close();
     fCntr++; 
  } 

}
