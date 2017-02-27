#ifndef Analysis_TriggerStudies_L1TriggersSemiLep_h
#define Analysis_TriggerStudies_L1TriggersSemiLep_h 1

#include "Analysis/Tools/interface/Analysis.h"

using namespace std;
using namespace analysis;
using namespace analysis::tools;


bool L1SingleMu3(Analysis & );
bool L1SingleJet20 (Analysis &);
bool L1DoubleJet32Eta2p2_Mu8_DEtaMax0p4_DPhiMax0p4_DEta1p6(Analysis &);
bool L1DoubleJet34Eta2p2_Mu8_DEtaMax0p4_DPhiMax0p4_DEta1p6(Analysis &);
bool L1DoubleJet36Eta2p2_Mu8_DEtaMax0p4_DPhiMax0p4_DEta1p6(Analysis &);
bool L1DoubleJet32Eta2p2_Mu10_DEtaMax0p4_DPhiMax0p4_DEta1p6(Analysis &);
bool L1DoubleJet34Eta2p2_Mu10_DEtaMax0p4_DPhiMax0p4_DEta1p6(Analysis &);
bool L1DoubleJet36Eta2p2_Mu10_DEtaMax0p4_DPhiMax0p4_DEta1p6(Analysis &);
bool L1DoubleJet32Eta2p2_Mu12_DEtaMax0p4_DPhiMax0p4_DEta1p6(Analysis &);
bool L1DoubleJet34Eta2p2_Mu12_DEtaMax0p4_DPhiMax0p4_DEta1p6(Analysis &);
bool L1DoubleJet36Eta2p2_Mu12_DEtaMax0p4_DPhiMax0p4_DEta1p6(Analysis &);

float detaL1Mumax = 0.4;
float dphiL1Mumax = 0.4;
float detaL1Jetmax = 1.6;

// L1 triggers

// ----------------------------------------------------------------------
bool L1SingleMu3(Analysis & analysis)
{
   if ( ! analysis.triggerResult("HLT_L1SingleMu3_v") ) return false;
   
   return true;
}

bool L1SingleJet20(Analysis & analysis)
{
   if ( ! analysis.triggerResult("HLT_L1SingleJet20_v") ) return false;
   
   return true;
}


//---------------------------------------------------------------------------------------- DOUBLEJET34 MU varies
bool L1DoubleJet34Eta2p2_Mu10_dEtaMuMax0p4_dPhiMuMax0p4_dEtaDoubleJet1p6(Analysis & analysis)
{
  if ( ! L1SingleMu3(analysis)   ) return false;
  if ( ! L1SingleJet20(analysis) ) return false;

  // Muon pt > 10 and eta < 2.2 selection at L1 
  auto l1Mu3s = analysis.collection<TriggerObject>("hltL1sSingleMu3");
  if ( l1Mu3s->size() < 1 )        return false; // additional consistency check
  
  std::vector<TriggerObject> selL1Mu10s;
  for ( int m = 0 ; m < l1Mu3s->size() ; ++m )
      {
         TriggerObject l1Mu3 = l1Mu3s->at(m);
         if ( l1Mu3.pt() >= 10 && fabs(l1Mu3.eta()) <= 2.2 ) 
         {
            selL1Mu10s.push_back(l1Mu3);
	    //   std::cout << l1Mu3.pt() << std::endl;
         }
         
      }
  if ( selL1Mu10s.size() < 1 )  return false;
  //) std::cout << "Mu10 Selection DONE " << std::endl;
   
   // Dijet selection with pt > 34 and eta < 2.2  
   auto l1Jet20s = analysis.collection<TriggerObject>("hltL1sSingleJet20");
   //  if ( l1Jet20s->size() < 4 )        return false; // since JET duplicated collection!!!!!!!!! DEFAULT 2
   cout << "Jetsize  " << l1Jet20s->size() << endl;
   std::vector<TriggerObject> selL1Jet30s;  
   for ( int j = 0 ; j < l1Jet20s->size() ; ++j )
      {
	TriggerObject l1Jet20 = l1Jet20s->at(j);
 
        if ( j != l1Jet20s->size()-1 ){   // last has no next
	TriggerObject nextJet = l1Jet20s->at(j+1);
	if ( (nextJet.pt() - l1Jet20.pt()) < 1e-7 && (nextJet.eta() - l1Jet20.eta()) < 1e-7 && (nextJet.phi() - l1Jet20.phi()) < 1e-7 ) continue ; // SKIP JET DUPLICATED
        }
       	cout <<  l1Jet20.pt() << " " << l1Jet20.eta()  << endl;
	if ( l1Jet20.pt() >= 34 && fabs(l1Jet20.eta()) <= 2.2 ) 
         {
            selL1Jet30s.push_back(l1Jet20);
         }
      } 
   if ( selL1Jet30s.size() < 2 || selL1Jet30s.size() > 5 )  return false;  //
   //   std::cout << "Dijet30 Selection DONE " << std::endl;

    // Muon matching to at least one Jet within [0.4,0,4] eta-phi square 
       int MuonMatches = 0;
       for ( int m = 0 ; m < abs(selL1Mu10s.size()) ; ++m )
         {
	    for ( int j = 0 ; j < abs(selL1Jet30s.size()) ; ++j )
	      { 
		float detaL1MuJet =  fabs(selL1Mu10s[m].eta() - selL1Jet30s[j].eta()) ;
		float dphiL1MuJet =  fabs(selL1Mu10s[m].phi() - selL1Jet30s[j].phi()) ;
		
		if ( detaL1MuJet < detaL1Mumax && dphiL1MuJet < dphiL1Mumax ) 
		  {
		    MuonMatches += 1;  //just need one muon match per jet
		    break;
		  }
	       }
          }
       if ( MuonMatches == 0 )      return false;
       //  std::cout << "# MuonMatches " << MuonMatches << std::endl;
       
     // At least a pair of well-separted jets with DEta < 1.6
       int JetDetaMatches = 0;
       for ( int j = 0 ; j <  abs(selL1Jet30s.size()) ; ++j )
         {
	    for ( int jj = 0 ; jj < abs(selL1Jet30s.size()) ; ++jj )
	      { 
		if (j == jj) continue;
		float detaL1JetJet =  fabs(selL1Jet30s[j].eta() - selL1Jet30s[jj].eta()) ;
		//	float dphiL1MuJet =  fabs(selL1Mu10s[m].phi() - selL1Jet30s[j].phi()) ;
		
		if ( detaL1JetJet < detaL1Jetmax ) 
		  {
		    JetDetaMatches += 1;  //just need one muon match
		    break;
		  }
	       }
          }
       if ( JetDetaMatches == 0 )    return false;
       //  std::cout << "# JetDetaMatches " << JetDetaMatches << std::endl; 
  
     return true;

     //can ADD ASYMMETRIC PT 

}

bool L1DoubleJet34Eta2p2_Mu12_dEtaMuMax0p4_dPhiMuMax0p4_dEtaDoubleJet1p6(Analysis & analysis)
{
  if ( ! L1SingleMu3(analysis)   ) return false;
  if ( ! L1SingleJet20(analysis) ) return false;

  // Muon pt > 10 and eta < 2.2 selection at L1 
  auto l1Mu3s = analysis.collection<TriggerObject>("hltL1sSingleMu3");
  if ( l1Mu3s->size() < 1 )        return false; // additional consistency check
  
  std::vector<TriggerObject> selL1Mu10s;
  for ( int m = 0 ; m < l1Mu3s->size() ; ++m )
      {
         TriggerObject l1Mu3 = l1Mu3s->at(m);
         if ( l1Mu3.pt() >= 12 && fabs(l1Mu3.eta()) <= 2.2 ) 
         {
            selL1Mu10s.push_back(l1Mu3);
    //   std::cout << l1Mu3.pt() << std::endl;
         }
         
      }
  if ( selL1Mu10s.size() < 1 )  return false;
  //) std::cout << "Mu10 Selection DONE " << std::endl;
   
   // Dijet selection with pt > 34 and eta < 2.2  
   auto l1Jet20s = analysis.collection<TriggerObject>("hltL1sSingleJet20");
   //  if ( l1Jet20s->size() < 4 )        return false; // since JET duplicated collection!!!!!!!!! DEFAULT 2

   std::vector<TriggerObject> selL1Jet30s;  
   for ( int j = 0 ; j < l1Jet20s->size() ; ++j )
      {
	TriggerObject l1Jet20 = l1Jet20s->at(j);
 
        if ( j != l1Jet20s->size()-1 ){   // last has no next
	TriggerObject nextJet = l1Jet20s->at(j+1);
	if ( (nextJet.pt() - l1Jet20.pt()) < 1e-7 && (nextJet.eta() - l1Jet20.eta()) < 1e-7 && (nextJet.phi() - l1Jet20.phi()) < 1e-7 ) continue ; // SKIP JET DUPLICATED
        }

	if ( l1Jet20.pt() >= 34 && fabs(l1Jet20.eta()) <= 2.2 ) 
         {
            selL1Jet30s.push_back(l1Jet20);
         }
      } 
   if ( selL1Jet30s.size() < 2 || selL1Jet30s.size() > 5 )  return false;  //
   //   std::cout << "Dijet30 Selection DONE " << std::endl;

    // Muon matching to at least one Jet within [0.4,0,4] eta-phi square 
       int MuonMatches = 0;
       for ( int m = 0 ; m < abs(selL1Mu10s.size()) ; ++m )
         {
	    for ( int j = 0 ; j < abs(selL1Jet30s.size()) ; ++j )
	      { 
		float detaL1MuJet =  fabs(selL1Mu10s[m].eta() - selL1Jet30s[j].eta()) ;
		float dphiL1MuJet =  fabs(selL1Mu10s[m].phi() - selL1Jet30s[j].phi()) ;
		
		if ( detaL1MuJet < detaL1Mumax && dphiL1MuJet < dphiL1Mumax ) 
		  {
		    MuonMatches += 1;  //just need one muon match per jet
		    break;
		  }
	       }
          }
       if ( MuonMatches == 0 )      return false;
       //  std::cout << "# MuonMatches " << MuonMatches << std::endl;
       
     // At least a pair of well-separted jets with DEta < 1.6
       int JetDetaMatches = 0;
       for ( int j = 0 ; j <  abs(selL1Jet30s.size()) ; ++j )
         {
	    for ( int jj = 0 ; jj < abs(selL1Jet30s.size()) ; ++jj )
	      { 
		if (j == jj) continue;
		float detaL1JetJet =  fabs(selL1Jet30s[j].eta() - selL1Jet30s[jj].eta()) ;
		//	float dphiL1MuJet =  fabs(selL1Mu10s[m].phi() - selL1Jet30s[j].phi()) ;
		
		if ( detaL1JetJet < detaL1Jetmax ) 
		  {
		    JetDetaMatches += 1;  //just need one muon match
		    break;
		  }
	       }
          }
       if ( JetDetaMatches == 0 )    return false;
       //  std::cout << "# JetDetaMatches " << JetDetaMatches << std::endl; 
  
     return true;

     //can ADD ASYMMETRIC PT 

}

bool L1DoubleJet34Eta2p2_Mu8_dEtaMuMax0p4_dPhiMuMax0p4_dEtaDoubleJet1p6(Analysis & analysis)
{
  if ( ! L1SingleMu3(analysis)   ) return false;
  if ( ! L1SingleJet20(analysis) ) return false;

  // Muon pt > 10 and eta < 2.2 selection at L1 
  auto l1Mu3s = analysis.collection<TriggerObject>("hltL1sSingleMu3");
  if ( l1Mu3s->size() < 1 )        return false; // additional consistency check
  
  std::vector<TriggerObject> selL1Mu10s;
  for ( int m = 0 ; m < l1Mu3s->size() ; ++m )
      {
         TriggerObject l1Mu3 = l1Mu3s->at(m);
         if ( l1Mu3.pt() >= 8 && fabs(l1Mu3.eta()) <= 2.2 ) 
         {
            selL1Mu10s.push_back(l1Mu3);
    //   std::cout << l1Mu3.pt() << std::endl;
         }
         
      }
  if ( selL1Mu10s.size() < 1 )  return false;
  //) std::cout << "Mu10 Selection DONE " << std::endl;
   
   // Dijet selection with pt > 34 and eta < 2.2  
   auto l1Jet20s = analysis.collection<TriggerObject>("hltL1sSingleJet20");
   //  if ( l1Jet20s->size() < 4 )        return false; // since JET duplicated collection!!!!!!!!! DEFAULT 2

   std::vector<TriggerObject> selL1Jet30s;  
   for ( int j = 0 ; j < l1Jet20s->size() ; ++j )
      {
	TriggerObject l1Jet20 = l1Jet20s->at(j);
 
        if ( j != l1Jet20s->size()-1 ){   // last has no next
	TriggerObject nextJet = l1Jet20s->at(j+1);
	if ( (nextJet.pt() - l1Jet20.pt()) < 1e-7 && (nextJet.eta() - l1Jet20.eta()) < 1e-7 && (nextJet.phi() - l1Jet20.phi()) < 1e-7 ) continue ; // SKIP JET DUPLICATED
        }

	if ( l1Jet20.pt() >= 34 && fabs(l1Jet20.eta()) <= 2.2 ) 
         {
            selL1Jet30s.push_back(l1Jet20);
         }
      } 
   if ( selL1Jet30s.size() < 2 || selL1Jet30s.size() > 5 )  return false;  //
   //   std::cout << "Dijet30 Selection DONE " << std::endl;

    // Muon matching to at least one Jet within [0.4,0,4] eta-phi square 
       int MuonMatches = 0;
       for ( int m = 0 ; m < abs(selL1Mu10s.size()) ; ++m )
         {
	    for ( int j = 0 ; j < abs(selL1Jet30s.size()) ; ++j )
	      { 
		float detaL1MuJet =  fabs(selL1Mu10s[m].eta() - selL1Jet30s[j].eta()) ;
		float dphiL1MuJet =  fabs(selL1Mu10s[m].phi() - selL1Jet30s[j].phi()) ;
		
		if ( detaL1MuJet < detaL1Mumax && dphiL1MuJet < dphiL1Mumax ) 
		  {
		    MuonMatches += 1;  //just need one muon match per jet
		    break;
		  }
	       }
          }
       if ( MuonMatches == 0 )      return false;
       //  std::cout << "# MuonMatches " << MuonMatches << std::endl;
       
     // At least a pair of well-separted jets with DEta < 1.6
       int JetDetaMatches = 0;
       for ( int j = 0 ; j <  abs(selL1Jet30s.size()) ; ++j )
         {
	    for ( int jj = 0 ; jj < abs(selL1Jet30s.size()) ; ++jj )
	      { 
		if (j == jj) continue;
		float detaL1JetJet =  fabs(selL1Jet30s[j].eta() - selL1Jet30s[jj].eta()) ;
		//	float dphiL1MuJet =  fabs(selL1Mu10s[m].phi() - selL1Jet30s[j].phi()) ;
		
		if ( detaL1JetJet < detaL1Jetmax ) 
		  {
		    JetDetaMatches += 1;  //just need one muon match
		    break;
		  }
	       }
          }
       if ( JetDetaMatches == 0 )    return false;
       //  std::cout << "# JetDetaMatches " << JetDetaMatches << std::endl; 
  
     return true;

     //can ADD ASYMMETRIC PT 

}

//---------------------------------------------------------------------------------------- DOUBLEJET32 MU varies





bool L1DoubleJet32Eta2p2_Mu10_dEtaMuMax0p4_dPhiMuMax0p4_dEtaDoubleJet1p6(Analysis & analysis)
{
  if ( ! L1SingleMu3(analysis)   ) return false;
  if ( ! L1SingleJet20(analysis) ) return false;

  // Muon pt > 10 and eta < 2.2 selection at L1 
  auto l1Mu3s = analysis.collection<TriggerObject>("hltL1sSingleMu3");
  if ( l1Mu3s->size() < 1 )        return false; // additional consistency check
  
  std::vector<TriggerObject> selL1Mu10s;
  for ( int m = 0 ; m < l1Mu3s->size() ; ++m )
      {
         TriggerObject l1Mu3 = l1Mu3s->at(m);
         if ( l1Mu3.pt() >= 10 && fabs(l1Mu3.eta()) <= 2.2 ) 
         {
            selL1Mu10s.push_back(l1Mu3);
    //   std::cout << l1Mu3.pt() << std::endl;
         }
         
      }
  if ( selL1Mu10s.size() < 1 )  return false;
  //) std::cout << "Mu10 Selection DONE " << std::endl;
   
   // Dijet selection with pt > 34 and eta < 2.2  
   auto l1Jet20s = analysis.collection<TriggerObject>("hltL1sSingleJet20");
   //  if ( l1Jet20s->size() < 4 )        return false; // since JET duplicated collection!!!!!!!!! DEFAULT 2

   std::vector<TriggerObject> selL1Jet30s;  
   for ( int j = 0 ; j < l1Jet20s->size() ; ++j )
      {
	TriggerObject l1Jet20 = l1Jet20s->at(j);
 
        if ( j != l1Jet20s->size()-1 ){   // last has no next
	TriggerObject nextJet = l1Jet20s->at(j+1);
	if ( (nextJet.pt() - l1Jet20.pt()) < 1e-7 && (nextJet.eta() - l1Jet20.eta()) < 1e-7 && (nextJet.phi() - l1Jet20.phi()) < 1e-7 ) continue ; // SKIP JET DUPLICATED
        }

	if ( l1Jet20.pt() >= 32 && fabs(l1Jet20.eta()) <= 2.2 ) 
         {
            selL1Jet30s.push_back(l1Jet20);
         }
      } 
   if ( selL1Jet30s.size() < 2 || selL1Jet30s.size() > 5 )  return false;  //
   //   std::cout << "Dijet30 Selection DONE " << std::endl;

    // Muon matching to at least one Jet within [0.4,0,4] eta-phi square 
       int MuonMatches = 0;
       for ( int m = 0 ; m < abs(selL1Mu10s.size()) ; ++m )
         {
	    for ( int j = 0 ; j < abs(selL1Jet30s.size()) ; ++j )
	      { 
		float detaL1MuJet =  fabs(selL1Mu10s[m].eta() - selL1Jet30s[j].eta()) ;
		float dphiL1MuJet =  fabs(selL1Mu10s[m].phi() - selL1Jet30s[j].phi()) ;
		
		if ( detaL1MuJet < detaL1Mumax && dphiL1MuJet < dphiL1Mumax ) 
		  {
		    MuonMatches += 1;  //just need one muon match per jet
		    break;
		  }
	       }
          }
       if ( MuonMatches == 0 )      return false;
       //  std::cout << "# MuonMatches " << MuonMatches << std::endl;
       
     // At least a pair of well-separted jets with DEta < 1.6
       int JetDetaMatches = 0;
       for ( int j = 0 ; j <  abs(selL1Jet30s.size()) ; ++j )
         {
	    for ( int jj = 0 ; jj < abs(selL1Jet30s.size()) ; ++jj )
	      { 
		if (j == jj) continue;
		float detaL1JetJet =  fabs(selL1Jet30s[j].eta() - selL1Jet30s[jj].eta()) ;
		//	float dphiL1MuJet =  fabs(selL1Mu10s[m].phi() - selL1Jet30s[j].phi()) ;
		
		if ( detaL1JetJet < detaL1Jetmax ) 
		  {
		    JetDetaMatches += 1;  //just need one muon match
		    break;
		  }
	       }
          }
       if ( JetDetaMatches == 0 )    return false;
       //  std::cout << "# JetDetaMatches " << JetDetaMatches << std::endl; 
  
     return true;

     //can ADD ASYMMETRIC PT 

}

bool L1DoubleJet32Eta2p2_Mu12_dEtaMuMax0p4_dPhiMuMax0p4_dEtaDoubleJet1p6(Analysis & analysis)
{
  if ( ! L1SingleMu3(analysis)   ) return false;
  if ( ! L1SingleJet20(analysis) ) return false;

  // Muon pt > 10 and eta < 2.2 selection at L1 
  auto l1Mu3s = analysis.collection<TriggerObject>("hltL1sSingleMu3");
  if ( l1Mu3s->size() < 1 )        return false; // additional consistency check
  
  std::vector<TriggerObject> selL1Mu10s;
  for ( int m = 0 ; m < l1Mu3s->size() ; ++m )
      {
         TriggerObject l1Mu3 = l1Mu3s->at(m);
         if ( l1Mu3.pt() >= 12 && fabs(l1Mu3.eta()) <= 2.2 ) 
         {
            selL1Mu10s.push_back(l1Mu3);
    //   std::cout << l1Mu3.pt() << std::endl;
         }
         
      }
  if ( selL1Mu10s.size() < 1 )  return false;
  //) std::cout << "Mu10 Selection DONE " << std::endl;
   
   // Dijet selection with pt > 34 and eta < 2.2  
   auto l1Jet20s = analysis.collection<TriggerObject>("hltL1sSingleJet20");
   //  if ( l1Jet20s->size() < 4 )        return false; // since JET duplicated collection!!!!!!!!! DEFAULT 2

   std::vector<TriggerObject> selL1Jet30s;  
   for ( int j = 0 ; j < l1Jet20s->size() ; ++j )
      {
	TriggerObject l1Jet20 = l1Jet20s->at(j);
 
        if ( j != l1Jet20s->size()-1 ){   // last has no next
	TriggerObject nextJet = l1Jet20s->at(j+1);
	if ( (nextJet.pt() - l1Jet20.pt()) < 1e-7 && (nextJet.eta() - l1Jet20.eta()) < 1e-7 && (nextJet.phi() - l1Jet20.phi()) < 1e-7 ) continue ; // SKIP JET DUPLICATED
        }

	if ( l1Jet20.pt() >= 32 && fabs(l1Jet20.eta()) <= 2.2 ) 
         {
            selL1Jet30s.push_back(l1Jet20);
         }
      } 
   if ( selL1Jet30s.size() < 2 || selL1Jet30s.size() > 5 )  return false;  //
   //   std::cout << "Dijet30 Selection DONE " << std::endl;

    // Muon matching to at least one Jet within [0.4,0,4] eta-phi square 
       int MuonMatches = 0;
       for ( int m = 0 ; m < abs(selL1Mu10s.size()) ; ++m )
         {
	    for ( int j = 0 ; j < abs(selL1Jet30s.size()) ; ++j )
	      { 
		float detaL1MuJet =  fabs(selL1Mu10s[m].eta() - selL1Jet30s[j].eta()) ;
		float dphiL1MuJet =  fabs(selL1Mu10s[m].phi() - selL1Jet30s[j].phi()) ;
		
		if ( detaL1MuJet < detaL1Mumax && dphiL1MuJet < dphiL1Mumax ) 
		  {
		    MuonMatches += 1;  //just need one muon match per jet
		    break;
		  }
	       }
          }
       if ( MuonMatches == 0 )      return false;
       //  std::cout << "# MuonMatches " << MuonMatches << std::endl;
       
     // At least a pair of well-separted jets with DEta < 1.6
       int JetDetaMatches = 0;
       for ( int j = 0 ; j <  abs(selL1Jet30s.size()) ; ++j )
         {
	    for ( int jj = 0 ; jj < abs(selL1Jet30s.size()) ; ++jj )
	      { 
		if (j == jj) continue;
		float detaL1JetJet =  fabs(selL1Jet30s[j].eta() - selL1Jet30s[jj].eta()) ;
		//	float dphiL1MuJet =  fabs(selL1Mu10s[m].phi() - selL1Jet30s[j].phi()) ;
		
		if ( detaL1JetJet < detaL1Jetmax ) 
		  {
		    JetDetaMatches += 1;  //just need one muon match
		    break;
		  }
	       }
          }
       if ( JetDetaMatches == 0 )    return false;
       //  std::cout << "# JetDetaMatches " << JetDetaMatches << std::endl; 
  
     return true;

     //can ADD ASYMMETRIC PT 

}

bool L1DoubleJet32Eta2p2_Mu8_dEtaMuMax0p4_dPhiMuMax0p4_dEtaDoubleJet1p6(Analysis & analysis)
{
  if ( ! L1SingleMu3(analysis)   ) return false;
  if ( ! L1SingleJet20(analysis) ) return false;

  // Muon pt > 10 and eta < 2.2 selection at L1 
  auto l1Mu3s = analysis.collection<TriggerObject>("hltL1sSingleMu3");
  if ( l1Mu3s->size() < 1 )        return false; // additional consistency check
  
  std::vector<TriggerObject> selL1Mu10s;
  for ( int m = 0 ; m < l1Mu3s->size() ; ++m )
      {
         TriggerObject l1Mu3 = l1Mu3s->at(m);
         if ( l1Mu3.pt() >= 8 && fabs(l1Mu3.eta()) <= 2.2 ) 
         {
            selL1Mu10s.push_back(l1Mu3);
    //   std::cout << l1Mu3.pt() << std::endl;
         }
         
      }
  if ( selL1Mu10s.size() < 1 )  return false;
  //) std::cout << "Mu10 Selection DONE " << std::endl;
   
   // Dijet selection with pt > 34 and eta < 2.2  
   auto l1Jet20s = analysis.collection<TriggerObject>("hltL1sSingleJet20");
   //  if ( l1Jet20s->size() < 4 )        return false; // since JET duplicated collection!!!!!!!!! DEFAULT 2

   std::vector<TriggerObject> selL1Jet30s;  
   for ( int j = 0 ; j < l1Jet20s->size() ; ++j )
      {
	TriggerObject l1Jet20 = l1Jet20s->at(j);
 
        if ( j != l1Jet20s->size()-1 ){   // last has no next
	TriggerObject nextJet = l1Jet20s->at(j+1);
	if ( (nextJet.pt() - l1Jet20.pt()) < 1e-7 && (nextJet.eta() - l1Jet20.eta()) < 1e-7 && (nextJet.phi() - l1Jet20.phi()) < 1e-7 ) continue ; // SKIP JET DUPLICATED
        }

	if ( l1Jet20.pt() >= 32 && fabs(l1Jet20.eta()) <= 2.2 ) 
         {
            selL1Jet30s.push_back(l1Jet20);
         }
      } 
   if ( selL1Jet30s.size() < 2 || selL1Jet30s.size() > 5 )  return false;  //
   //   std::cout << "Dijet30 Selection DONE " << std::endl;

    // Muon matching to at least one Jet within [0.4,0,4] eta-phi square 
       int MuonMatches = 0;
       for ( int m = 0 ; m < abs(selL1Mu10s.size()) ; ++m )
         {
	    for ( int j = 0 ; j < abs(selL1Jet30s.size()) ; ++j )
	      { 
		float detaL1MuJet =  fabs(selL1Mu10s[m].eta() - selL1Jet30s[j].eta()) ;
		float dphiL1MuJet =  fabs(selL1Mu10s[m].phi() - selL1Jet30s[j].phi()) ;
		
		if ( detaL1MuJet < detaL1Mumax && dphiL1MuJet < dphiL1Mumax ) 
		  {
		    MuonMatches += 1;  //just need one muon match per jet
		    break;
		  }
	       }
          }
       if ( MuonMatches == 0 )      return false;
       //  std::cout << "# MuonMatches " << MuonMatches << std::endl;
       
     // At least a pair of well-separted jets with DEta < 1.6
       int JetDetaMatches = 0;
       for ( int j = 0 ; j <  abs(selL1Jet30s.size()) ; ++j )
         {
	    for ( int jj = 0 ; jj < abs(selL1Jet30s.size()) ; ++jj )
	      { 
		if (j == jj) continue;
		float detaL1JetJet =  fabs(selL1Jet30s[j].eta() - selL1Jet30s[jj].eta()) ;
		//	float dphiL1MuJet =  fabs(selL1Mu10s[m].phi() - selL1Jet30s[j].phi()) ;
		
		if ( detaL1JetJet < detaL1Jetmax ) 
		  {
		    JetDetaMatches += 1;  //just need one muon match
		    break;
		  }
	       }
          }
       if ( JetDetaMatches == 0 )    return false;
       //  std::cout << "# JetDetaMatches " << JetDetaMatches << std::endl; 
  
     return true;

     //can ADD ASYMMETRIC PT 

}

//---------------------------------------------------------------------------------------- DOUBLEJET36 MU varies



bool L1DoubleJet36Eta2p2_Mu10_dEtaMuMax0p4_dPhiMuMax0p4_dEtaDoubleJet1p6(Analysis & analysis)
{
  if ( ! L1SingleMu3(analysis)   ) return false;
  if ( ! L1SingleJet20(analysis) ) return false;

  // Muon pt > 10 and eta < 2.2 selection at L1 
  auto l1Mu3s = analysis.collection<TriggerObject>("hltL1sSingleMu3");
  if ( l1Mu3s->size() < 1 )        return false; // additional consistency check
  
  std::vector<TriggerObject> selL1Mu10s;
  for ( int m = 0 ; m < l1Mu3s->size() ; ++m )
      {
         TriggerObject l1Mu3 = l1Mu3s->at(m);
         if ( l1Mu3.pt() >= 10 && fabs(l1Mu3.eta()) <= 2.2 ) 
         {
            selL1Mu10s.push_back(l1Mu3);
    //   std::cout << l1Mu3.pt() << std::endl;
         }
         
      }
  if ( selL1Mu10s.size() < 1 )  return false;
  //) std::cout << "Mu10 Selection DONE " << std::endl;
   
   // Dijet selection with pt > 34 and eta < 2.2  
   auto l1Jet20s = analysis.collection<TriggerObject>("hltL1sSingleJet20");
   //  if ( l1Jet20s->size() < 4 )        return false; // since JET duplicated collection!!!!!!!!! DEFAULT 2

   std::vector<TriggerObject> selL1Jet30s;  
   for ( int j = 0 ; j < l1Jet20s->size() ; ++j )
      {
	TriggerObject l1Jet20 = l1Jet20s->at(j);
 
        if ( j != l1Jet20s->size()-1 ){   // last has no next
	TriggerObject nextJet = l1Jet20s->at(j+1);
	if ( (nextJet.pt() - l1Jet20.pt()) < 1e-7 && (nextJet.eta() - l1Jet20.eta()) < 1e-7 && (nextJet.phi() - l1Jet20.phi()) < 1e-7 ) continue ; // SKIP JET DUPLICATED
        }

	if ( l1Jet20.pt() >= 36 && fabs(l1Jet20.eta()) <= 2.2 ) 
         {
            selL1Jet30s.push_back(l1Jet20);
         }
      } 
   if ( selL1Jet30s.size() < 2 || selL1Jet30s.size() > 5 )  return false;  //
   //   std::cout << "Dijet30 Selection DONE " << std::endl;

    // Muon matching to at least one Jet within [0.4,0,4] eta-phi square 
       int MuonMatches = 0;
       for ( int m = 0 ; m < abs(selL1Mu10s.size()) ; ++m )
         {
	    for ( int j = 0 ; j < abs(selL1Jet30s.size()) ; ++j )
	      { 
		float detaL1MuJet =  fabs(selL1Mu10s[m].eta() - selL1Jet30s[j].eta()) ;
		float dphiL1MuJet =  fabs(selL1Mu10s[m].phi() - selL1Jet30s[j].phi()) ;
		
		if ( detaL1MuJet < detaL1Mumax && dphiL1MuJet < dphiL1Mumax ) 
		  {
		    MuonMatches += 1;  //just need one muon match per jet
		    break;
		  }
	       }
          }
       if ( MuonMatches == 0 )      return false;
       //  std::cout << "# MuonMatches " << MuonMatches << std::endl;
       
     // At least a pair of well-separted jets with DEta < 1.6
       int JetDetaMatches = 0;
       for ( int j = 0 ; j <  abs(selL1Jet30s.size()) ; ++j )
         {
	    for ( int jj = 0 ; jj < abs(selL1Jet30s.size()) ; ++jj )
	      { 
		if (j == jj) continue;
		float detaL1JetJet =  fabs(selL1Jet30s[j].eta() - selL1Jet30s[jj].eta()) ;
		//	float dphiL1MuJet =  fabs(selL1Mu10s[m].phi() - selL1Jet30s[j].phi()) ;
		
		if ( detaL1JetJet < detaL1Jetmax ) 
		  {
		    JetDetaMatches += 1;  //just need one muon match
		    break;
		  }
	       }
          }
       if ( JetDetaMatches == 0 )    return false;
       //  std::cout << "# JetDetaMatches " << JetDetaMatches << std::endl; 
  
     return true;

     //can ADD ASYMMETRIC PT 

}

bool L1DoubleJet36Eta2p2_Mu12_dEtaMuMax0p4_dPhiMuMax0p4_dEtaDoubleJet1p6(Analysis & analysis)
{
  if ( ! L1SingleMu3(analysis)   ) return false;
  if ( ! L1SingleJet20(analysis) ) return false;

  // Muon pt > 10 and eta < 2.2 selection at L1 
  auto l1Mu3s = analysis.collection<TriggerObject>("hltL1sSingleMu3");
  if ( l1Mu3s->size() < 1 )        return false; // additional consistency check
  
  std::vector<TriggerObject> selL1Mu10s;
  for ( int m = 0 ; m < l1Mu3s->size() ; ++m )
      {
         TriggerObject l1Mu3 = l1Mu3s->at(m);
         if ( l1Mu3.pt() >= 12 && fabs(l1Mu3.eta()) <= 2.2 ) 
         {
            selL1Mu10s.push_back(l1Mu3);
    //   std::cout << l1Mu3.pt() << std::endl;
         }
         
      }
  if ( selL1Mu10s.size() < 1 )  return false;
  //) std::cout << "Mu10 Selection DONE " << std::endl;
   
   // Dijet selection with pt > 34 and eta < 2.2  
   auto l1Jet20s = analysis.collection<TriggerObject>("hltL1sSingleJet20");
   //  if ( l1Jet20s->size() < 4 )        return false; // since JET duplicated collection!!!!!!!!! DEFAULT 2

   std::vector<TriggerObject> selL1Jet30s;  
   for ( int j = 0 ; j < l1Jet20s->size() ; ++j )
      {
	TriggerObject l1Jet20 = l1Jet20s->at(j);
 
        if ( j != l1Jet20s->size()-1 ){   // last has no next
	TriggerObject nextJet = l1Jet20s->at(j+1);
	if ( (nextJet.pt() - l1Jet20.pt()) < 1e-7 && (nextJet.eta() - l1Jet20.eta()) < 1e-7 && (nextJet.phi() - l1Jet20.phi()) < 1e-7 ) continue ; // SKIP JET DUPLICATED
        }

	if ( l1Jet20.pt() >= 36 && fabs(l1Jet20.eta()) <= 2.2 ) 
         {
            selL1Jet30s.push_back(l1Jet20);
         }
      } 
   if ( selL1Jet30s.size() < 2 || selL1Jet30s.size() > 5 )  return false;  //
   //   std::cout << "Dijet30 Selection DONE " << std::endl;

    // Muon matching to at least one Jet within [0.4,0,4] eta-phi square 
       int MuonMatches = 0;
       for ( int m = 0 ; m < abs(selL1Mu10s.size()) ; ++m )
         {
	    for ( int j = 0 ; j < abs(selL1Jet30s.size()) ; ++j )
	      { 
		float detaL1MuJet =  fabs(selL1Mu10s[m].eta() - selL1Jet30s[j].eta()) ;
		float dphiL1MuJet =  fabs(selL1Mu10s[m].phi() - selL1Jet30s[j].phi()) ;
		
		if ( detaL1MuJet < detaL1Mumax && dphiL1MuJet < dphiL1Mumax ) 
		  {
		    MuonMatches += 1;  //just need one muon match per jet
		    break;
		  }
	       }
          }
       if ( MuonMatches == 0 )      return false;
       //  std::cout << "# MuonMatches " << MuonMatches << std::endl;
       
     // At least a pair of well-separted jets with DEta < 1.6
       int JetDetaMatches = 0;
       for ( int j = 0 ; j <  abs(selL1Jet30s.size()) ; ++j )
         {
	    for ( int jj = 0 ; jj < abs(selL1Jet30s.size()) ; ++jj )
	      { 
		if (j == jj) continue;
		float detaL1JetJet =  fabs(selL1Jet30s[j].eta() - selL1Jet30s[jj].eta()) ;
		//	float dphiL1MuJet =  fabs(selL1Mu10s[m].phi() - selL1Jet30s[j].phi()) ;
		
		if ( detaL1JetJet < detaL1Jetmax ) 
		  {
		    JetDetaMatches += 1;  //just need one muon match
		    break;
		  }
	       }
          }
       if ( JetDetaMatches == 0 )    return false;
       //  std::cout << "# JetDetaMatches " << JetDetaMatches << std::endl; 
  
     return true;

     //can ADD ASYMMETRIC PT 

}

bool L1DoubleJet36Eta2p2_Mu8_dEtaMuMax0p4_dPhiMuMax0p4_dEtaDoubleJet1p6(Analysis & analysis)
{
  if ( ! L1SingleMu3(analysis)   ) return false;
  if ( ! L1SingleJet20(analysis) ) return false;

  // Muon pt > 10 and eta < 2.2 selection at L1 
  auto l1Mu3s = analysis.collection<TriggerObject>("hltL1sSingleMu3");
  if ( l1Mu3s->size() < 1 )        return false; // additional consistency check
  
  std::vector<TriggerObject> selL1Mu10s;
  for ( int m = 0 ; m < l1Mu3s->size() ; ++m )
      {
         TriggerObject l1Mu3 = l1Mu3s->at(m);
         if ( l1Mu3.pt() >= 8 && fabs(l1Mu3.eta()) <= 2.2 ) 
         {
            selL1Mu10s.push_back(l1Mu3);
    //   std::cout << l1Mu3.pt() << std::endl;
         }
         
      }
  if ( selL1Mu10s.size() < 1 )  return false;
  //) std::cout << "Mu10 Selection DONE " << std::endl;
   
   // Dijet selection with pt > 34 and eta < 2.2  
   auto l1Jet20s = analysis.collection<TriggerObject>("hltL1sSingleJet20");
   //  if ( l1Jet20s->size() < 4 )        return false; // since JET duplicated collection!!!!!!!!! DEFAULT 2

   std::vector<TriggerObject> selL1Jet30s;  
   for ( int j = 0 ; j < l1Jet20s->size() ; ++j )
      {
	TriggerObject l1Jet20 = l1Jet20s->at(j);
 
        if ( j != l1Jet20s->size()-1 ){   // last has no next
	TriggerObject nextJet = l1Jet20s->at(j+1);
	if ( (nextJet.pt() - l1Jet20.pt()) < 1e-7 && (nextJet.eta() - l1Jet20.eta()) < 1e-7 && (nextJet.phi() - l1Jet20.phi()) < 1e-7 ) continue ; // SKIP JET DUPLICATED
        }

	if ( l1Jet20.pt() >= 36 && fabs(l1Jet20.eta()) <= 2.2 ) 
         {
            selL1Jet30s.push_back(l1Jet20);
         }
      } 
   if ( selL1Jet30s.size() < 2 || selL1Jet30s.size() > 5 )  return false;  //
   //   std::cout << "Dijet30 Selection DONE " << std::endl;

    // Muon matching to at least one Jet within [0.4,0,4] eta-phi square 
       int MuonMatches = 0;
       for ( int m = 0 ; m < abs(selL1Mu10s.size()) ; ++m )
         {
	    for ( int j = 0 ; j < abs(selL1Jet30s.size()) ; ++j )
	      { 
		float detaL1MuJet =  fabs(selL1Mu10s[m].eta() - selL1Jet30s[j].eta()) ;
		float dphiL1MuJet =  fabs(selL1Mu10s[m].phi() - selL1Jet30s[j].phi()) ;
		
		if ( detaL1MuJet < detaL1Mumax && dphiL1MuJet < dphiL1Mumax ) 
		  {
		    MuonMatches += 1;  //just need one muon match per jet
		    break;
		  }
	       }
          }
       if ( MuonMatches == 0 )      return false;
       //  std::cout << "# MuonMatches " << MuonMatches << std::endl;
       
     // At least a pair of well-separted jets with DEta < 1.6
       int JetDetaMatches = 0;
       for ( int j = 0 ; j <  abs(selL1Jet30s.size()) ; ++j )
         {
	    for ( int jj = 0 ; jj < abs(selL1Jet30s.size()) ; ++jj )
	      { 
		if (j == jj) continue;
		float detaL1JetJet =  fabs(selL1Jet30s[j].eta() - selL1Jet30s[jj].eta()) ;
		//	float dphiL1MuJet =  fabs(selL1Mu10s[m].phi() - selL1Jet30s[j].phi()) ;
		
		if ( detaL1JetJet < detaL1Jetmax ) 
		  {
		    JetDetaMatches += 1;  //just need one muon match
		    break;
		  }
	       }
          }
       if ( JetDetaMatches == 0 )    return false;
       //  std::cout << "# JetDetaMatches " << JetDetaMatches << std::endl; 
  
     return true;

     //can ADD ASYMMETRIC PT 

}



//DOUBLE COUNTING OF JETS, MUON ARE FINE

#endif  // Analysis_TriggerStudies_L1TriggersSemiLep_h
