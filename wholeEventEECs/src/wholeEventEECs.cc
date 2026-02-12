#include "wholeEventEECs.h"

#include <format>


wholeEventEECs::wholeEventEECs(const std::string &name)
  : SubsysReco(name)
{
}

bool wholeEventEECs::isGoodTrigger(PHCompositeNode *topNode)
{
    bool goodTrigger = false;
    TriggerAnalyzer *ta = new TriggerAnalyzer();
    ta->UseEmulator(false);
    ta->decodeTriggers(topNode);
    goodTrigger = ta->didTriggerFire("Jet 10 GeV");
    return goodTrigger;
}

bool wholeEventEECs::isGoodDijet()
{
    if(vtx)
    {
        m_vtx_z = vtx->get_z();
        if(std::abs(m_vtx_z) > 60.0)
        {
            std::cout << "   vertex not in range" << std::endl;
            return false;
        };
    }
    else
    {
        m_vtx_z = 0.0;
    }

    if(jets->size() < 2)
    {
        std::cout << "   fewer than 2 jets" << std::endl;
        return false;
    }

    Jet *leadJet = nullptr;
    Jet *subleadJet = nullptr;
    float lead_pT = 0.0;
    float sublead_pT = 0.0;
    
    for(int i=0; i<(int)jets->size(); i++)
    {
        Jet *j = jets->get_jet(i);
        float pT = j->get_pt();
        if(pT > lead_pT)
        {
            if(leadJet)
            {
                subleadJet = leadJet;
                sublead_pT = lead_pT;
            }
            leadJet = j;
            lead_pT = pT;
        }
    }
    
    if(!leadJet || !subleadJet)
    {
        std::cout << "   bad lead or sublead jet" << std::endl;
        return false;
    }

    if(lead_pT < sublead_pT)
    {
        std::cout << "   lead pT < sub pT" << std::endl;
        return false;
    }

    if(lead_pT < 20.9 || sublead_pT < 9.4)
    {
        std::cout << "   lead or sub pT too low: lead < 20.9? " << lead_pT << "   sub < 9.4? " << sublead_pT << std::endl;
        return false;
    }

    if(std::abs(leadJet->get_eta()) > 0.7 || std::abs(subleadJet->get_eta()) > 0.7)
    {
        std::cout << "   lead or sub not in right eta range" << std::endl;
        return false;
    }

    float dPhi = subleadJet->get_phi() - leadJet->get_phi();
    if(dPhi > M_PI) dPhi -= 2*M_PI;
    if(dPhi < -M_PI) dPhi += 2*M_PI;
    if(std::abs(dPhi) < 0.75*M_PI)
    {
        std::cout << "   dPhi between sub and lead too small (should be greater than " << 0.75*M_PI << "): " << std::abs(dPhi) << std::endl;
        return false;
    }

    return true;
}

std::array<float, 3> wholeEventEECs::correct_for_vertex(std::array<float, 3> center)
{
    std::array<float, 3> correctedTower;

    float z = center[2]*sinh(center[0]);
    if(m_vtx_z != 0.0)
    {
        z += m_vtx_z;
        float newEta = asinh(z/center[2]);
        if(!std::isnan(newEta) && !std::isinf(newEta))
        {
            correctedTower = {newEta, center[1], center[2]};
        }
    }
    else
    {
        correctedTower = {center[0], center[1], center[2]};
    }

    return correctedTower;

}

std::array<float, 5> wholeEventEECs::calc_observables(std::map<std::array<float, 3>, float>::iterator &iterA, std::map<std::array<float, 3>, float>::iterator &iterB)
{
    float dEta = std::abs(iterA->first[0] - iterB->first[0]);
    float dPhi = iterA->first[1] - iterB->first[1];
    if(dPhi > M_PI) dPhi += 2*M_PI;
    if(dPhi < -M_PI) dPhi -= 2*M_PI;
    dPhi = std::abs(dPhi);

    float dR = sqrt(dEta*dEta + dPhi*dPhi);


    float pTA = iterA->second/cosh(iterA->first[0]);
    float pxA = pTA*cos(iterA->first[1]);
    float pyA = pTA*sin(iterA->first[1]);
    float pzA = pTA*sinh(iterA->first[0]);

    float pTB = iterB->second/cosh(iterB->first[0]);
    float pxB = pTB*cos(iterB->first[1]);
    float pyB = pTB*sin(iterB->first[1]);
    float pzB = pTB*sinh(iterB->first[0]);

    float cosTheta = (pxA*pxB + pyA*pyB + pzA*pzB)/(iterA->second * iterB->second);

    return std::array<float, 5> {dR, dEta, dPhi, std::acos(cosTheta), (((float)1.0)-cosTheta)/((float)2.0)};
}


int wholeEventEECs::InitRun(PHCompositeNode *topNode)
{

    emcalTowerInfoContainer = findNode::getClass<TowerInfoContainer>(topNode,"TOWERINFO_CALIB_CEMC_RETOWER");
    ihcalTowerInfoContainer = findNode::getClass<TowerInfoContainer>(topNode,"TOWERINFO_CALIB_HCALIN");
    ohcalTowerInfoContainer = findNode::getClass<TowerInfoContainer>(topNode,"TOWERINFO_CALIB_HCALOUT");

    if(!emcalTowerInfoContainer || !ihcalTowerInfoContainer || !ohcalTowerInfoContainer)
    {
        std::cout << "One or more TowerInfoContainers missing. Exiting" << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
    }

    emcal_geom = findNode::getClass<RawTowerGeomContainer_Cylinderv1>(topNode, "TOWERGEOM_CEMC");
    ihcal_geom = findNode::getClass<RawTowerGeomContainer_Cylinderv1>(topNode, "TOWERGEOM_HCALIN");
    ohcal_geom = findNode::getClass<RawTowerGeomContainer_Cylinderv1>(topNode, "TOWERGEOM_HCALOUT");
    if(!emcal_geom || !ihcal_geom || !ohcal_geom)
    {
        std::cout << "One or more Tower Geometry Containers missing. Exiting" << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
    }

    vtxMap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
    if(!vtxMap)
    {
        vtxMap = findNode::getClass<GlobalVertexMap>(topNode, "MbdVertexMap");
        if(!vtxMap)
        {
            std::cout << "No global vertex map found. Exiting" << std::endl;
            return Fun4AllReturnCodes::ABORTRUN;
        }
    }

    //jets = findNode::getClass<JetContainer>(topNode, "AntiKt_TowerInfo_r04");
    jets = findNode::getClass<JetContainer>(topNode, "AntiKt_unsubtracted_r04");
    if(!jets)
    {
        std::cout << "AntiKt_TowerInfo_r04 jet node missing. Exiting" << std::endl;
        return Fun4AllReturnCodes::ABORTRUN;
    }

    m_outputFile = new TFile(m_outputName.c_str(),"RECREATE");

    nEvents = new TH1D("nEvents","Number of dijet events",1,0.5,1.5);

    for(int i=0; i<5; i++)
    {
        for(int j=0; j<4; j++)
        {
            EEC[i][j] = new TH1D(std::format("EEC_{}_{}", obsNames[i], caloNames[j]).c_str(), std::format("Whole Event EEC {};{}", caloNames[j], obsSymbols[i]).c_str(),50,0.0,obsMax[i]);
            EEC_corr[i][j] = new TH1D(std::format("EEC_corr_{}_{}", obsNames[i], caloNames[j]).c_str(), std::format("Vertex Corrected Whole Event EEC {};{}", caloNames[j], obsSymbols[i]).c_str(),50,0.0,obsMax[i]);
        }

        /*
        EEC_RL[i] = new TH1D(std::format("EEC_RL_{}", caloNames[i]).c_str(), std::format("Whole Event EEC {};#DeltaR", caloNames[i]).c_str(),50,0.0,sqrt(2.2*2.2 + M_PI*M_PI));
        EEC_eta[i] = new TH1D(std::format("EEC_eta_{}", caloNames[i]).c_str(), std::format("Whole Event EEC {};#Delta#eta", caloNames[i]).c_str(),50,0.0,2.2);
        EEC_phi[i] = new TH1D(std::format("EEC_phi_{}", caloNames[i]).c_str(), std::format("Whole Event EEC {};#Delta#phi", caloNames[i]).c_str(),50,0.0,M_PI);
        EEC_theta[i] = new TH1D(std::format("EEC_theta_{}", caloNames[i]).c_str(), std::format("Whole Event EEC {};#theta", caloNames[i]).c_str(),50,0.0,M_PI);
        EEC_z[i] = new TH1D(std::format("EEC_z_{}", caloNames[i]).c_str(), std::format("Whole Event EEC {};z=(1-cos#theta)/2", caloNames[i]).c_str(),50,0.0,1.0);
    
        EEC_RL_corr[i] = new TH1D(std::format("EEC_RL_corr_{}", caloNames[i]).c_str(), std::format("Vertex Corrected Whole Event EEC {};#DeltaR", caloNames[i]).c_str(),50,0.0,sqrt(2.2*2.2 + M_PI*M_PI));
        EEC_eta_corr[i] = new TH1D(std::format("EEC_eta_corr_{}", caloNames[i]).c_str(), std::format("Vertex Corrected Whole Event EEC {};#Delta#eta", caloNames[i]).c_str(),50,0.0,2.2);
        EEC_phi_corr[i] = new TH1D(std::format("EEC_phi_corr_{}", caloNames[i]).c_str(), std::format("Vertex Corrected Whole Event EEC {};#Delta#phi", caloNames[i]).c_str(),50,0.0,M_PI);
        EEC_theta_corr[i] = new TH1D(std::format("EEC_theta_corr_{}", caloNames[i]).c_str(), std::format("Vertex Corrected Whole Event EEC {};#theta", caloNames[i]).c_str(),50,0.0,M_PI);
        EEC_z_corr[i] = new TH1D(std::format("EEC_z_corr_{}", caloNames[i]).c_str(), std::format("Vertex Corrected Whole Event EEC {};z=(1-cos#theta)/2", caloNames[i]).c_str(),50,0.0,1.0);
        */
    }


    return Fun4AllReturnCodes::EVENT_OK;

}

int wholeEventEECs::process_event(PHCompositeNode *topNode)
{

    EMCalTowers.clear();
    IHCalTowers.clear();
    OHCalTowers.clear();
    AllTowers.clear();

    EMCalCorrectedTowers.clear();
    IHCalCorrectedTowers.clear();
    OHCalCorrectedTowers.clear();
    AllCorrectedTowers.clear();

    m_eventIndex++;


    std::cout << "working on event " << m_eventIndex << ". There have been " << m_nGoodEvents << " good events so far." << std::endl;

    bool goodTrigger = isGoodTrigger(topNode);
    if(!goodTrigger)
    {
        std::cout << "   trigger not found" << std::endl;
        return Fun4AllReturnCodes::EVENT_OK;
    }

    if(vtxMap->empty())
    {
        std::cout << "      empty vertex map, setting vz to origin" << std::endl;
        vtx = nullptr;
    }
    else
    {
        for(auto vItr : *vtxMap)
        {
            if(vItr.first == 0 || vItr.first == GlobalVertex::VTXTYPE::MBD || vItr.first == GlobalVertex::VTXTYPE::SVTX || vItr.first == GlobalVertex::VTXTYPE::SVTX_MBD)
            {
                vtx = vItr.second;
            }
        }
    }

    bool goodDijet = isGoodDijet();
    if(!goodDijet)
    {
        std::cout << "   good dijet not found" << std::endl;
        return Fun4AllReturnCodes::EVENT_OK;
    }

    m_nGoodEvents++;
    nEvents->Fill(1);

    for(int i=0; i<(int) emcalTowerInfoContainer->size(); i++)
    {
        if(!emcalTowerInfoContainer->get_tower_at_channel(i)->get_isGood()) continue;
        emcal_geom->set_calorimeter_id(RawTowerDefs::CEMC);
        ihcal_geom->set_calorimeter_id(RawTowerDefs::HCALIN);
        auto key = emcalTowerInfoContainer->encode_key(i);
        auto tower = emcalTowerInfoContainer->get_tower_at_channel(i);
        int phiBin = emcalTowerInfoContainer->getTowerPhiBin(key);
        int etaBin = emcalTowerInfoContainer->getTowerEtaBin(key);
        float phiCenter = ihcal_geom->get_phicenter(phiBin);
        float etaCenter = ihcal_geom->get_etacenter(etaBin);
        float r = emcal_geom->get_radius();
        std::array<float, 3> center {etaCenter, phiCenter, r};
        EMCalTowers.insert(std::make_pair(center, tower->get_energy()));
        if(AllTowers.find(center) != AllTowers.end()) AllTowers.at(center) += tower->get_energy();
        else AllTowers[center] = tower->get_energy();

        std::array<float, 3> newCenter = correct_for_vertex(center);
        EMCalCorrectedTowers.insert(std::make_pair(newCenter, tower->get_energy()));
        if(AllCorrectedTowers.find(newCenter) != AllCorrectedTowers.end()) AllCorrectedTowers.at(newCenter) += tower->get_energy();
        else AllCorrectedTowers[newCenter] = tower->get_energy();
    }

    for(int i=0; i<(int) ihcalTowerInfoContainer->size(); i++)
    {
        if(!ihcalTowerInfoContainer->get_tower_at_channel(i)->get_isGood()) continue;
        ihcal_geom->set_calorimeter_id(RawTowerDefs::HCALIN);
        auto key = ihcalTowerInfoContainer->encode_key(i);
        auto tower = ihcalTowerInfoContainer->get_tower_at_channel(i);
        int phiBin = ihcalTowerInfoContainer->getTowerPhiBin(key);
        int etaBin = ihcalTowerInfoContainer->getTowerEtaBin(key);
        float phiCenter = ihcal_geom->get_phicenter(phiBin);
        float etaCenter = ihcal_geom->get_etacenter(etaBin);
        float r = ihcal_geom->get_radius();
        std::array<float, 3> center {etaCenter, phiCenter, r};
        IHCalTowers.insert(std::make_pair(center, tower->get_energy()));
        if(AllTowers.find(center) != AllTowers.end()) AllTowers.at(center) += tower->get_energy();
        else AllTowers[center] = tower->get_energy();

        std::array<float, 3> newCenter = correct_for_vertex(center);
        IHCalCorrectedTowers.insert(std::make_pair(newCenter, tower->get_energy()));
        if(AllCorrectedTowers.find(newCenter) != AllCorrectedTowers.end()) AllCorrectedTowers.at(newCenter) += tower->get_energy();
        else AllCorrectedTowers[newCenter] = tower->get_energy();
    }

    for(int i=0; i<(int) ohcalTowerInfoContainer->size(); i++)
    {
        if(!ohcalTowerInfoContainer->get_tower_at_channel(i)->get_isGood()) continue;
        ohcal_geom->set_calorimeter_id(RawTowerDefs::HCALOUT);
        auto key = ohcalTowerInfoContainer->encode_key(i);
        auto tower = ohcalTowerInfoContainer->get_tower_at_channel(i);
        int phiBin = ohcalTowerInfoContainer->getTowerPhiBin(key);
        int etaBin = ohcalTowerInfoContainer->getTowerEtaBin(key);
        float phiCenter = ohcal_geom->get_phicenter(phiBin);
        float etaCenter = ohcal_geom->get_etacenter(etaBin);
        float r = ohcal_geom->get_radius();
        std::array<float, 3> center {etaCenter, phiCenter, r};
        OHCalTowers.insert(std::make_pair(center, tower->get_energy()));
        if(AllTowers.find(center) != AllTowers.end()) AllTowers.at(center) += tower->get_energy();
        else AllTowers[center] = tower->get_energy();

        std::array<float, 3> newCenter = correct_for_vertex(center);
        OHCalCorrectedTowers.insert(std::make_pair(newCenter, tower->get_energy()));
        if(AllCorrectedTowers.find(newCenter) != AllCorrectedTowers.end()) AllCorrectedTowers.at(newCenter) += tower->get_energy();
        else AllCorrectedTowers[newCenter] = tower->get_energy();
    }

    //no correcting for vertex position
    for(auto it1 = EMCalTowers.begin(); it1 != EMCalTowers.end(); ++it1)
    {
        float pT1 = it1->second / cosh(it1->first[0]);
        if(pT1 < 0.3) continue;
        for(auto it2 = std::next(it1); it2 != EMCalTowers.end(); ++it2)
        {
            float pT2 = it2->second / cosh(it2->first[0]);
            if(pT2 < 0.3) continue;


            std::array<float, 5> observables = calc_observables(it1, it2);
            for(int i=0; i<5; i++)
            {
                EEC[i][0]->Fill(observables[i], pT1 * pT2);
            }
        }
    }

    for(auto it1 = IHCalTowers.begin(); it1 != IHCalTowers.end(); ++it1)
    {
        float pT1 = it1->second / cosh(it1->first[0]);
        if(pT1 < 0.3) continue; 
        for(auto it2 = std::next(it1); it2 != IHCalTowers.end(); ++it2)
        {
            float pT2 = it2->second / cosh(it2->first[0]);
            if(pT2 < 0.3) continue;


            std::array<float, 5> observables = calc_observables(it1, it2);
            for(int i=0; i<5; i++)
            {
                EEC[i][1]->Fill(observables[i], pT1 * pT2);
            }
        }
    }

    for(auto it1 = OHCalTowers.begin(); it1 != OHCalTowers.end(); ++it1)
    {
        float pT1 = it1->second / cosh(it1->first[0]);
        if(pT1 < 0.3) continue;
        for(auto it2 = std::next(it1); it2 != OHCalTowers.end(); ++it2)
        {
            float pT2 = it2->second / cosh(it2->first[0]);
            if(pT2 < 0.3) continue;


            std::array<float, 5> observables = calc_observables(it1, it2);
            for(int i=0; i<5; i++)
            {
                EEC[i][2]->Fill(observables[i], pT1 * pT2);
            }
        }
    }

    for(auto it1 = AllTowers.begin(); it1 != AllTowers.end(); ++it1)
    {
        float pT1 = it1->second / cosh(it1->first[0]);
        if(pT1 < 0.3) continue;
        for(auto it2 = std::next(it1); it2 != AllTowers.end(); ++it2)
        {
            float pT2 = it2->second / cosh(it2->first[0]);
            if(pT2 < 0.3) continue;


            std::array<float, 5> observables = calc_observables(it1, it2);
            for(int i=0; i<5; i++)
            {
                EEC[i][3]->Fill(observables[i], pT1 * pT2);
            }
        }
    }

    //correcting for vertex position
    for(auto it1 = EMCalCorrectedTowers.begin(); it1 != EMCalCorrectedTowers.end(); ++it1)
    {
        float pT1 = it1->second / cosh(it1->first[0]);
        if(pT1 < 0.3) continue;
        for(auto it2 = std::next(it1); it2 != EMCalCorrectedTowers.end(); ++it2)
        {
            float pT2 = it2->second / cosh(it2->first[0]);
            if(pT2 < 0.3) continue;


            std::array<float, 5> observables = calc_observables(it1, it2);
            for(int i=0; i<5; i++)
            {
                EEC_corr[i][0]->Fill(observables[i], pT1 * pT2);
            }
        }
    }

    for(auto it1 = IHCalCorrectedTowers.begin(); it1 != IHCalCorrectedTowers.end(); ++it1)
    {
        float pT1 = it1->second / cosh(it1->first[0]);
        if(pT1 < 0.3) continue;
        for(auto it2 = std::next(it1); it2 != IHCalCorrectedTowers.end(); ++it2)
        {
            float pT2 = it2->second / cosh(it2->first[0]);
            if(pT2 < 0.3) continue;


            std::array<float, 5> observables = calc_observables(it1, it2);
            for(int i=0; i<5; i++)
            {
                EEC_corr[i][1]->Fill(observables[i], pT1 * pT2);
            }
        }
    }

    for(auto it1 = OHCalCorrectedTowers.begin(); it1 != OHCalCorrectedTowers.end(); ++it1)
    {
        float pT1 = it1->second / cosh(it1->first[0]);
        if(pT1 < 0.3) continue;
        for(auto it2 = std::next(it1); it2 != OHCalCorrectedTowers.end(); ++it2)
        {
            float pT2 = it2->second / cosh(it2->first[0]);
            if(pT2 < 0.3) continue;


            std::array<float, 5> observables = calc_observables(it1, it2);
            for(int i=0; i<5; i++)
            {
                EEC_corr[i][2]->Fill(observables[i], pT1 * pT2);
            }
        }
    }

    for(auto it1 = AllCorrectedTowers.begin(); it1 != AllCorrectedTowers.end(); ++it1)
    {
        float pT1 = it1->second / cosh(it1->first[0]);
        if(pT1 < 0.3) continue;
        for(auto it2 = std::next(it1); it2 != AllCorrectedTowers.end(); ++it2)
        {
            float pT2 = it2->second / cosh(it2->first[0]);
            if(pT2 < 0.3) continue;


            std::array<float, 5> observables = calc_observables(it1, it2);
            for(int i=0; i<5; i++)
            {
                EEC_corr[i][3]->Fill(observables[i], pT1 * pT2);
            }
        }
    }

    return Fun4AllReturnCodes::EVENT_OK;
}

int wholeEventEECs::End(PHCompositeNode */*topNode*/)
{
    m_outputFile->cd();
    nEvents->Write();
    for(int i=0; i<5; i++)
    {
        for(int j=0; j<4; j++)
        {
            EEC[i][j]->Write();
            EEC_corr[i][j]->Write();
        }
    }
    m_outputFile->Close();

    return Fun4AllReturnCodes::EVENT_OK;
}
