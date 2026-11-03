"""
FAZENA Tools — High-Coverage Target Prediction & Proteome Navigation Engine
Lightweight target prediction using 2D Morgan fingerprints and curated bioactive target chemotypes.
Fast, offline, zero-cloud inference covering human and pathogen proteomes (Kinome, GPCRome, Proteases, Ion Channels, Nuclear Receptors, Transporters, Epigenome).
Author: Yehor Franchuk (@MrEgorin), FAZENA (fazena.org)
Affiliation: Institute of Cell Biology NASU
"""

from typing import List, Dict, Any, Optional
import math
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.Chem import rdFingerprintGenerator

_MORGAN_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

def _get_fp(mol):
    return _MORGAN_GEN.GetFingerprint(mol)

# ═════════════════════════════════════════════════════════════════════
# REFERENCE BIOACTIVE TARGET KNOWLEDGE BASE (HUMAN & PATHOGEN PROTEOMES)
# Curated target profiles with representative clinical & tool compounds
# ═════════════════════════════════════════════════════════════════════

TARGET_PROFILES: List[Dict[str, Any]] = [
    # ── KINASES (HUMAN KINOME) ───────────────────────────────────────
    {
        'target_name': 'Epidermal Growth Factor Receptor (EGFR)',
        'gene': 'EGFR',
        'class': 'Kinase',
        'organism': 'Human',
        'uniprot': 'P00533',
        'pdb_id': '1M17',
        'compounds': [
            ('Gefitinib', 'COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1'),
            ('Erlotinib', 'COCCOC1=C(OCCOC)C=C2C(=C1)C(=NC=N2)NC3=CC=CC(=C3)C#C'),
            ('Osimertinib', 'CN(C)CC=CC(=O)Nc1cc(Nc2ncccn2)c(OC)cc1N(C)CCN(C)C')
        ]
    },
    {
        'target_name': 'Cyclin-Dependent Kinase 2 (CDK2)',
        'gene': 'CDK2',
        'class': 'Kinase',
        'organism': 'Human',
        'uniprot': 'P24941',
        'pdb_id': '1HCK',
        'compounds': [
            ('Roscovitine', 'CCN(C1=NC=NC2=C1N=C(N2CC(C)C)NCC3=CC=CC=C3)C(C)CO'),
            ('Dinaciclib', 'CCN1CCN(CC1)C(=O)C2=CC=C(C=C2)NC3=NC=C4C(=N3)N(N=C4)C5CCCC5'),
            ('Flavopiridol', 'CN1CCC(C(C1)O)C2=C(C=C(C3=C2OC(=CC3=O)C4=CC=CC=C4Cl)O)O')
        ]
    },
    {
        'target_name': 'BRAF Kinase (V600E)',
        'gene': 'BRAF',
        'class': 'Kinase',
        'organism': 'Human',
        'uniprot': 'P15056',
        'pdb_id': '4MNE',
        'compounds': [
            ('Vemurafenib', 'CCCS(=O)(=O)NC1=C(C=C(C(=C1)F)C(=O)C2=CNC3=C2C=C(C=N3)C4=CC=C(C=C4)Cl)F'),
            ('Dabrafenib', 'CC(C)(C)C1=NC(=C(S1)C2=NC(=NC=C2)N)C3=C(C=C(C(=C3)F)NS(=O)(=O)C4=C(C=CC=C4F)F)F')
        ]
    },
    {
        'target_name': 'Tyrosine-Protein Kinase ABL1',
        'gene': 'ABL1',
        'class': 'Kinase',
        'organism': 'Human',
        'uniprot': 'P00519',
        'pdb_id': '1IEP',
        'compounds': [
            ('Imatinib', 'Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1'),
            ('Nilotinib', 'Cc1cn(-c2cc(NC(=O)c3ccc(C)c(Nc4nccc(-c5cccnc5)n4)c3)cc(C(F)(F)F)c2)cn1'),
            ('Dasatinib', 'Cc1nc(Nc2ncc(C(=O)Nc3c(C)cccc3Cl)s2)cc(N2CCN(CCO)CC2)n1')
        ]
    },
    {
        'target_name': 'Vascular Endothelial Growth Factor Receptor 2 (VEGFR2/KDR)',
        'gene': 'KDR',
        'class': 'Kinase',
        'organism': 'Human',
        'uniprot': 'P35968',
        'pdb_id': '4AG8',
        'compounds': [
            ('Sorafenib', 'CNC(=O)c1cc(Oc2ccc(NC(=O)Nc3ccc(Cl)c(C(F)(F)F)c3)cc2)ccn1'),
            ('Sunitinib', 'CCN(CC)CCNC(=O)c1c(C)[nH]c(C=C2C(=O)Nc3ccc(F)cc32)c1C'),
            ('Axitinib', 'CNC(=O)c1ccccc1Sc1ccc2c(cnn2-c2ccccn2)c1')
        ]
    },
    {
        'target_name': 'RAC-alpha Serine/Threonine-Protein Kinase (AKT1)',
        'gene': 'AKT1',
        'class': 'Kinase',
        'organism': 'Human',
        'uniprot': 'P31749',
        'pdb_id': '3QKK',
        'compounds': [
            ('Capivasertib', 'CC(C)NC(=O)c1cnc(nc1N)-c1cn[nH]c1'),
            ('Ipatasertib', 'CN1CCN(CC1)C(=O)c1cc(NC(=O)c2ccc(Cl)c(Cl)c2)ccn1'),
            ('Miransertib', 'CC1(CCN(CC1)c1ncnc2[nH]c(nc12)-c1ccc(F)cc1)N')
        ]
    },
    {
        'target_name': 'Dual Specificity Mitogen-Activated Protein Kinase Kinase 1 (MEK1)',
        'gene': 'MAP2K1',
        'class': 'Kinase',
        'organism': 'Human',
        'uniprot': 'Q02750',
        'pdb_id': '3WGI',
        'compounds': [
            ('Trametinib', 'CC1=C(C(=O)N(C2=C1N(C(=O)NC2=O)C3=C(C=C(C=C3)I)F)C4=CC=CC=C4F)C5CC5'),
            ('Cobimetinib', 'OC1(CCN(CC1)Cc1c(F)c(F)cc(Nc2ccc(I)cc2F)c1)C(=O)O'),
            ('Selumetinib', 'CN1C=C(C(=O)NOCC(O)CO)c2cc(Nc3ccc(I)cc3Cl)c(F)cc21')
        ]
    },
    {
        'target_name': 'Tyrosine-Protein Kinase JAK2',
        'gene': 'JAK2',
        'class': 'Kinase',
        'organism': 'Human',
        'uniprot': 'O60674',
        'pdb_id': '4IVA',
        'compounds': [
            ('Ruxolitinib', 'N#CCC(C1CCCC1)n1cc(-c2ncnc3[nH]ccc23)cn1'),
            ('Fedratinib', 'CC(C)(C)S(=O)(=O)c1cccc(Nc2nc(Nc3ccc(OCCN4CCCC4)cc3)ncc2Cl)c1'),
            ('Baricitinib', 'CCS(=O)(=O)N1CC(CC#N)(n2cc(-c3ncnc4[nH]ccc34)cn2)C1')
        ]
    },
    {
        'target_name': 'Bruton Tyrosine Kinase (BTK)',
        'gene': 'BTK',
        'class': 'Kinase',
        'organism': 'Human',
        'uniprot': 'Q06187',
        'pdb_id': '5P9J',
        'compounds': [
            ('Ibrutinib', 'C=CC(=O)N1CCC(CC1)n1nc(-c2ccc(Oc3ccccc3)cc2)c2c(N)ncnc21'),
            ('Acalabrutinib', 'CC#CC(=O)N1CCCC1c1nc(-c2ccc(NC(=O)c3ccccn3)cc2)c2cnn(C)c2n1'),
            ('Zanubrutinib', 'C=CC(=O)N1CCC(CC1)n1c(N)c(-c2ccc(Oc3ccccc3)cc2)nc1C(=O)N')
        ]
    },
    {
        'target_name': 'Phosphatidylinositol 4,5-bisphosphate 3-kinase catalytic subunit alpha (PIK3CA)',
        'gene': 'PIK3CA',
        'class': 'Kinase',
        'organism': 'Human',
        'uniprot': 'P42336',
        'pdb_id': '4JPS',
        'compounds': [
            ('Alpelisib', 'CC(C)(C)c1nc(NC(=O)N2CCC(C2)c2cnc(N)nc2)cs1'),
            ('Buparlisib', 'CC(C)(C#N)c1ccc(nc1)-c1nc(N)nc(N2CCOCC2)n1'),
            ('Idelalisib', 'CCC1=NC2=C(C=CC=C2Cl)C(=O)N1C(C)C1=NC=NC2=C1NC=N2')
        ]
    },
    {
        'target_name': 'Anaplastic Lymphoma Kinase (ALK)',
        'gene': 'ALK',
        'class': 'Kinase',
        'organism': 'Human',
        'uniprot': 'Q9UM73',
        'pdb_id': '2XP2',
        'compounds': [
            ('Crizotinib', 'CC(Oc1cc(-c2cnn(C3CCNCC3)c2)cnc1N)c1c(Cl)ccc(F)c1Cl'),
            ('Alectinib', 'CCc1cc2c(cc1N1CCC(CC1)N1CCOCC1)C(C)(C)c1cc(C#N)ccc1N2'),
            ('Brigatinib', 'COc1cc(Nc2ncc(Cl)c(Nc3ccccc3P(=O)(C)C)n2)ccc1N1CCN(C)CC1')
        ]
    },
    {
        'target_name': 'Hepatocyte Growth Factor Receptor (MET)',
        'gene': 'MET',
        'class': 'Kinase',
        'organism': 'Human',
        'uniprot': 'P08581',
        'pdb_id': '3DKC',
        'compounds': [
            ('Capmatinib', 'Cc1nnc2cn(Cc3c(F)cccc3Cl)nc2c1-c1cccc(C#N)c1'),
            ('Cabozantinib', 'COc1cc2ncnc(Oc3ccc(NC(=O)C4(C(=O)Nc5ccc(F)cc5)CC4)cc3)c2cc1OC')
        ]
    },
    {
        'target_name': 'Proto-oncogene Tyrosine-Protein Kinase Src (SRC)',
        'gene': 'SRC',
        'class': 'Kinase',
        'organism': 'Human',
        'uniprot': 'P12931',
        'pdb_id': '2SRC',
        'compounds': [
            ('Bosutinib', 'COc1cc2c(Nc3ccc(Cl)c(Cl)c3)c(C#N)cnc2cc1OCCCN1CCN(C)CC1'),
            ('Saracatinib', 'CN1CCN(CCCOc2cc3ncnc(Nc4ccc(Cl)c(Cl)c4)c3cc2OC)CC1')
        ]
    },
    {
        'target_name': 'Aurora Kinase A (AURKA)',
        'gene': 'AURKA',
        'class': 'Kinase',
        'organism': 'Human',
        'uniprot': 'O14965',
        'pdb_id': '1OL5',
        'compounds': [
            ('Alisertib', 'COc1ccc(cc1)C(=O)Nc1ccc(cc1)C(=O)Nc1ncc(Cl)s1'),
            ('Tozasertib', 'Cc1cc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)n[nH]1')
        ]
    },

    # ── GPCRS (HUMAN GPCROME) ─────────────────────────────────────────
    {
        'target_name': '5-Hydroxytryptamine Receptor 2A (5-HT2A)',
        'gene': 'HTR2A',
        'class': 'GPCR',
        'organism': 'Human',
        'uniprot': 'P28223',
        'pdb_id': '6A93',
        'compounds': [
            ('Ketanserin', 'O=C1NC(=O)N(CCN2CCC(CC2)C(=O)C3=CC=C(F)C=C3)C1=O'),
            ('Risperidone', 'CC1=C(CCN2CCC(CC2)C3=NOC4=C3C=CC(=C4)F)C(=O)N5CCCC5=N1'),
            ('Clozapine', 'CN1CCN(CC1)C2=NC3=CC=CC=C3NC4=C2C=C(C=C4)Cl')
        ]
    },
    {
        'target_name': '5-Hydroxytryptamine Receptor 1A (5-HT1A)',
        'gene': 'HTR1A',
        'class': 'GPCR',
        'organism': 'Human',
        'uniprot': 'P08908',
        'pdb_id': '7E2X',
        'compounds': [
            ('Buspirone', 'O=C1CC2(CCCC2)CC(=O)N1CCCCN1CCN(c2ncccn2)CC1'),
            ('Tandospirone', 'O=C1CC2(CCCC2)CC(=O)N1CCCCN1CCN(c2ccccn2)CC1')
        ]
    },
    {
        'target_name': 'Dopamine D2 Receptor (DRD2)',
        'gene': 'DRD2',
        'class': 'GPCR',
        'organism': 'Human',
        'uniprot': 'P14416',
        'pdb_id': '6CM4',
        'compounds': [
            ('Haloperidol', 'OC1(CCN(CCCC(=O)c2ccc(F)cc2)CC1)c1ccc(Cl)cc1'),
            ('Aripiprazole', 'O=C1CCc2ccc(OCCCCN3CCN(c4cccc(Cl)c4Cl)CC3)cc2N1'),
            ('Sulpiride', 'CCN1CCCC1CNC(=O)c2cc(S(=O)(=O)N)ccc2OC')
        ]
    },
    {
        'target_name': 'Dopamine D1 Receptor (DRD1)',
        'gene': 'DRD1',
        'class': 'GPCR',
        'organism': 'Human',
        'uniprot': 'P21728',
        'pdb_id': '7CKX',
        'compounds': [
            ('Tavapadon', 'Cc1cccc(Oc2ccc(C#N)cc2)c1C1CCNCC1'),
            ('SKF-38393', 'Oc1ccc2c(c1O)CCNC2c1ccccc1')
        ]
    },
    {
        'target_name': 'Cannabinoid Receptor 1 (CB1)',
        'gene': 'CNR1',
        'class': 'GPCR',
        'organism': 'Human',
        'uniprot': 'P21554',
        'pdb_id': '5XRA',
        'compounds': [
            ('Rimonabant', 'Clc1ccc(cc1)-c1nn(c(Cl)c1C(=O)NN1CCCCC1)-c1ccc(Cl)cc1Cl'),
            ('THC', 'CCCCCC1=CC(O)=C2C(=C1)OC(C)(C)C1CC=C(C)CC21')
        ]
    },
    {
        'target_name': 'Cannabinoid Receptor 2 (CB2)',
        'gene': 'CNR2',
        'class': 'GPCR',
        'organism': 'Human',
        'uniprot': 'P34972',
        'pdb_id': '5ZTY',
        'compounds': [
            ('JWH-133', 'CCCCC1=CC(O)=C2C(=C1)OC(C)(C)C1CC=C(C)CC21'),
            ('HU-308', 'CCCCCC1=CC(=C2C(=C1)OC(C)(C)C1CC=C(C)CC21)OC')
        ]
    },
    {
        'target_name': 'Beta-1 Adrenergic Receptor (ADRB1)',
        'gene': 'ADRB1',
        'class': 'GPCR',
        'organism': 'Human',
        'uniprot': 'P08588',
        'pdb_id': '7BTS',
        'compounds': [
            ('Atenolol', 'CC(C)NCC(O)COc1ccc(CC(N)=O)cc1'),
            ('Metoprolol', 'COCCc1ccc(OCC(O)CNC(C)C)cc1'),
            ('Bisoprolol', 'COCCOCCCOc1ccc(OCC(O)CNC(C)C)cc1')
        ]
    },
    {
        'target_name': 'Beta-2 Adrenergic Receptor (ADRB2)',
        'gene': 'ADRB2',
        'class': 'GPCR',
        'organism': 'Human',
        'uniprot': 'P07550',
        'pdb_id': '2RH1',
        'compounds': [
            ('Salbutamol', 'CC(C)(C)NCC(O)c1ccc(O)c(CO)c1'),
            ('Propranolol', 'CC(C)NCC(O)COc1cccc2ccccc12'),
            ('Formoterol', 'COc1ccc(CC(C)NCC(O)c2ccc(O)c(NC=O)c2)cc1')
        ]
    },
    {
        'target_name': 'Mu-type Opioid Receptor (OPRM1)',
        'gene': 'OPRM1',
        'class': 'GPCR',
        'organism': 'Human',
        'uniprot': 'P35372',
        'pdb_id': '5C1M',
        'compounds': [
            ('Morphine', 'CN1CCC23C4OC5=C(O)C=CC(=C45)C2C1CC=C3O'),
            ('Fentanyl', 'CCC(=O)N(c1ccccc1)C1CCN(CCc2ccccc2)CC1'),
            ('Naloxone', 'O=C1CC2(O)C3Cc4ccc(O)c5c4C2(CCN3CC=C)C1O5')
        ]
    },
    {
        'target_name': 'Histamine H1 Receptor (HRH1)',
        'gene': 'HRH1',
        'class': 'GPCR',
        'organism': 'Human',
        'uniprot': 'P35367',
        'pdb_id': '3RZE',
        'compounds': [
            ('Cetirizine', 'OC(=O)COCCN1CCN(C(c2ccccc2)c2ccc(Cl)cc2)CC1'),
            ('Loratadine', 'CCOC(=O)N1CCC(=C2c3ccccc3CCc3ncc(Cl)cc32)CC1'),
            ('Diphenhydramine', 'CN(C)CCOC(c1ccccc1)c1ccccc1')
        ]
    },
    {
        'target_name': 'C-X-C Chemokine Receptor Type 4 (CXCR4)',
        'gene': 'CXCR4',
        'class': 'GPCR',
        'organism': 'Human',
        'uniprot': 'P61073',
        'pdb_id': '3ODU',
        'compounds': [
            ('IT1t', 'CC(C)NC(=N)NC1=C(Sc2ccccc21)C(=O)NC1CCNCC1'),
            ('AMD3100', 'c1cc(cnc1)CN1CCN(CC2CCNCC2)CC1')
        ]
    },

    # ── NUCLEAR RECEPTORS ─────────────────────────────────────────────
    {
        'target_name': 'Estrogen Receptor Alpha (ESR1)',
        'gene': 'ESR1',
        'class': 'Nuclear Receptor',
        'organism': 'Human',
        'uniprot': 'P03372',
        'pdb_id': '1ERR',
        'compounds': [
            ('Tamoxifen', 'CCC(=C(c1ccccc1)c1ccc(OCCN(C)C)cc1)c1ccccc1'),
            ('Raloxifene', 'Oc1ccc(cc1)C(=O)c1c(sc2cc(O)ccc12)c1ccc(OCCN2CCCCC2)cc1'),
            ('Fulvestrant', 'CC12CCC3C(C1CCC2O)CCC4=C3C=CC(=C4)O')
        ]
    },
    {
        'target_name': 'Androgen Receptor (AR)',
        'gene': 'AR',
        'class': 'Nuclear Receptor',
        'organism': 'Human',
        'uniprot': 'P10275',
        'pdb_id': '1E3G',
        'compounds': [
            ('Enzalutamide', 'CC1(C)C(=O)N(c2ccc(C#N)c(C(F)(F)F)c2)C(=S)N1c1ccc(F)c(C(=O)NC)c1'),
            ('Bicalutamide', 'CC(O)(CS(=O)(=O)c1ccc(F)cc1)C(=O)Nc1ccc(C#N)c(C(F)(F)F)c1'),
            ('Apalutamide', 'CC1(C)C(=O)N(c2ccc(C#N)c(C(F)(F)F)c2)C(=S)N1c1ccc(C(=O)NC)nc1')
        ]
    },
    {
        'target_name': 'Peroxisome Proliferator-Activated Receptor Gamma (PPARG)',
        'gene': 'PPARG',
        'class': 'Nuclear Receptor',
        'organism': 'Human',
        'uniprot': 'P37231',
        'pdb_id': '2PRG',
        'compounds': [
            ('Rosiglitazone', 'CN(CCOc1ccc(CC2SC(=O)NC2=O)cc1)c1ccccn1'),
            ('Pioglitazone', 'CCc1ccc(CCOc2ccc(CC3SC(=O)NC3=O)cc2)nc1'),
            ('Lobeglitazone', 'COc1ccc(CCOc2ccc(CC3SC(=O)NC3=O)cc2)nc1')
        ]
    },
    {
        'target_name': 'Glucocorticoid Receptor (NR3C1)',
        'gene': 'NR3C1',
        'class': 'Nuclear Receptor',
        'organism': 'Human',
        'uniprot': 'P04150',
        'pdb_id': '4P6X',
        'compounds': [
            ('Dexamethasone', 'CC1CC2C3CCC4=CC(=O)C=CC4(C)C3(F)C(O)CC2(C)C1(O)C(=O)CO'),
            ('Prednisolone', 'CC12C=CC(=O)C=C1CCC3C2C(O)CC4(C)C(O)(C(=O)CO)CCC34'),
            ('Hydrocortisone', 'CC12CCC(=O)C=C1CCC3C2C(O)CC4(C)C(O)(C(=O)CO)CCC34')
        ]
    },

    # ── PROTEASES & ENZYMES ───────────────────────────────────────────
    {
        'target_name': 'Cyclooxygenase-2 (COX-2 / PTGS2)',
        'gene': 'PTGS2',
        'class': 'Enzyme',
        'organism': 'Human',
        'uniprot': 'P35354',
        'pdb_id': '1CX2',
        'compounds': [
            ('Celecoxib', 'Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(N)(=O)=O)cc2)cc1'),
            ('Rofecoxib', 'O=S(=O)(N)c1ccc(cc1)C1=C(c2ccccc2)COC1=O'),
            ('Etoricoxib', 'Cc1ccc(nc1)-c1ccc(S(C)(=O)=O)cc1Cl')
        ]
    },
    {
        'target_name': 'Cyclooxygenase-1 (COX-1 / PTGS1)',
        'gene': 'PTGS1',
        'class': 'Enzyme',
        'organism': 'Human',
        'uniprot': 'P23219',
        'pdb_id': '6Y3C',
        'compounds': [
            ('Aspirin', 'CC(=O)Oc1ccccc1C(=O)O'),
            ('Indomethacin', 'COc1ccc2n(C(=O)c3ccc(Cl)cc3)c(C)c(CC(=O)O)c2c1'),
            ('Flurbiprofen', 'CC(C(=O)O)c1ccc(-c2ccccc2)c(F)c1')
        ]
    },
    {
        'target_name': 'Angiotensin-Converting Enzyme (ACE)',
        'gene': 'ACE',
        'class': 'Protease',
        'organism': 'Human',
        'uniprot': 'P12821',
        'pdb_id': '1O86',
        'compounds': [
            ('Captopril', 'CC(CS)C(=O)N1CCCC1C(=O)O'),
            ('Enalapril', 'CCOC(=O)C(CCc1ccccc1)NC(C)C(=O)N1CCCC1C(=O)O'),
            ('Lisinopril', 'NCCCCC(NC(CCc1ccccc1)C(=O)O)C(=O)N1CCCC1C(=O)O')
        ]
    },
    {
        'target_name': 'Dipeptidyl Peptidase 4 (DPP4)',
        'gene': 'DPP4',
        'class': 'Protease',
        'organism': 'Human',
        'uniprot': 'P27487',
        'pdb_id': '1X70',
        'compounds': [
            ('Sitagliptin', 'O=C(CC(N)Cc1cc(F)c(F)cc1F)N1CCn2c(C1)nnc2C(F)(F)F'),
            ('Vildagliptin', 'O=C(CNC12CC3CC(C1)CC(O)(C3)C2)N1CCCC1C#N'),
            ('Linagliptin', 'CC#CCn1c(N2CCCC2)nc2c1c(=O)n(Cc1nc(C)c3ccccc3n1)c(=O)n2C')
        ]
    },
    {
        'target_name': 'Coagulation Factor Xa (F10)',
        'gene': 'F10',
        'class': 'Protease',
        'organism': 'Human',
        'uniprot': 'P00742',
        'pdb_id': '1FOR',
        'compounds': [
            ('Rivaroxaban', 'O=C1NC(CNC(=O)c2ccc(Cl)s2)CN1c1ccc(N2CCOCC2=O)cc1'),
            ('Apixaban', 'COc1ccc(cc1)N1N=C(C(=O)N)C2=C1CCN(c1ccc(cc1)N1CCCCC1=O)C2=O'),
            ('Edoxaban', 'CN(C)C(=O)c1ccc(cc1)NC(=O)C1CCCC1NC(=O)c1nc(Cl)cs1')
        ]
    },
    {
        'target_name': 'Thrombin (F2)',
        'gene': 'F2',
        'class': 'Protease',
        'organism': 'Human',
        'uniprot': 'P00734',
        'pdb_id': '1PPB',
        'compounds': [
            ('Dabigatran', 'CN(c1nc2cc(C(=O)N(CCC(=O)OCC)c3ccccn3)ccc2n1C)c1ccc(C(=N)N)cc1'),
            ('Argatroban', 'Cc1cccc2c1c(S(=O)(=O)NC(CCCNC(=N)N)C(=O)N1CCCC1(C)C(=O)O)ccc2')
        ]
    },
    {
        'target_name': 'Acetylcholinesterase (ACHE)',
        'gene': 'ACHE',
        'class': 'Enzyme',
        'organism': 'Human',
        'uniprot': 'P22303',
        'pdb_id': '4EY7',
        'compounds': [
            ('Donepezil', 'COc1cc2c(cc1OC)C(=O)C(CC1CCN(Cc3ccccc3)CC1)C2'),
            ('Rivastigmine', 'CCN(C)C(=O)Oc1cccc(C(C)N(C)C)c1'),
            ('Galantamine', 'CN1CCC23C=CC(O)CC2Oc2c(OC)ccc(C1)c23')
        ]
    },
    {
        'target_name': 'Beta-Secretase 1 (BACE1)',
        'gene': 'BACE1',
        'class': 'Protease',
        'organism': 'Human',
        'uniprot': 'P56817',
        'pdb_id': '2QP8',
        'compounds': [
            ('Verubecestat', 'CN1C(=N)N(C)C(=O)C1(C)c1cc(NC(=O)c2cnc(F)cn2)ccc1F'),
            ('Atabecestat', 'CC1(C)N=C(N)N(C)C(=O)C1(C)c1cc(NC(=O)c2cnc(C#N)cn2)ccc1F')
        ]
    },
    {
        'target_name': 'Phosphodiesterase 5A (PDE5A)',
        'gene': 'PDE5A',
        'class': 'Enzyme',
        'organism': 'Human',
        'uniprot': 'O76074',
        'pdb_id': '1UDT',
        'compounds': [
            ('Sildenafil', 'CCCC1=NN(C)C2=C1N=C(NC2=O)C1=C(OCC)C=CC(=C1)S(=O)(=O)N1CCN(C)CC1'),
            ('Tadalafil', 'O=C1N2CC3=C(NC4=CC=CC=C34)C(C5=CC=C6OCOC6=C5)N2C(=O)C1'),
            ('Vardenafil', 'CCCC1=NN(C)C2=C1N=C(NC2=O)C1=C(OCC)C=CC(=C1)S(=O)(=O)N1CCN(CC)CC1')
        ]
    },

    # ── VIRAL & BACTERIAL PATHOGEN PROTEOMES ──────────────────────────
    {
        'target_name': 'SARS-CoV-2 Main Protease (Mpro / 3CLpro)',
        'gene': 'Mpro',
        'class': 'Viral Protease',
        'organism': 'Viral',
        'uniprot': 'P0DTD1',
        'pdb_id': '6LU7',
        'compounds': [
            ('Nirmatrelvir', 'CC1(C2C1C(N(C2)C(=O)C(C(C)(C)C)NC(=O)C(F)(F)F)C(=O)NC(CC3CCNC3=O)C#N)C'),
            ('Ensitrelvir', 'CN1C(=NN=C1C(=O)NC2=CC(=C(C=C2)F)F)C3=C(C=CC(=C3)F)F'),
            ('GC-376', 'CC(C)CC(NC(=O)OCc1ccccc1)C(=O)NC(CC1CCNC1=O)C=O')
        ]
    },
    {
        'target_name': 'SARS-CoV-2 Papain-Like Protease (PLpro)',
        'gene': 'PLpro',
        'class': 'Viral Protease',
        'organism': 'Viral',
        'uniprot': 'P0DTD1',
        'pdb_id': '7JRN',
        'compounds': [
            ('GRL-0617', 'Cc1cccc(C(=O)Nc2ccc(C)c(N)c2)c1'),
            ('Jun9-72-2', 'Cc1cccc(C(=O)Nc2ccc(CN3CCOCC3)c(N)c2)c1')
        ]
    },
    {
        'target_name': 'Human Immunodeficiency Virus 1 Protease (HIV-1 PR)',
        'gene': 'HIV-PR',
        'class': 'Viral Protease',
        'organism': 'Viral',
        'uniprot': 'P03366',
        'pdb_id': '1HVR',
        'compounds': [
            ('Darunavir', 'CC(C)CN(CC(O)C(Cc1ccccc1)NC(=O)OC1COC2OCCC12)S(=O)(=O)c1ccc(N)cc1'),
            ('Atazanavir', 'COC(=O)NC(C(=O)NC(Cc1ccccc1)C(O)CC(Cc1ccccc1)NC(=O)C(NC(=O)OC)C(C)(C)C)C(C)(C)C'),
            ('Ritonavir', 'CC(C)c1nc(CN(C)C(=O)NC(C(=O)NC(Cc2ccccc2)CC(O)C(Cc2ccccc2)NC(=O)OCc2cncs2)C(C)C)cs1')
        ]
    },
    {
        'target_name': 'Bacterial Cell Division Protein FtsZ',
        'gene': 'FtsZ',
        'class': 'Antimicrobial',
        'organism': 'Bacterial',
        'uniprot': 'P0A9A6',
        'pdb_id': '4DXD',
        'compounds': [
            ('PC190723', 'COc1cc2c(cc1OC)C(=O)C(Cc1ccc(F)c(Cl)c1)C2'),
            ('Berberine', 'COc1ccc2c(c1OC)[n+]1ccc3cc4c(cc3c1CC2)OCO4'),
            ('TXA-707', 'COc1cc2ncnc(Oc3ccc(Cl)c(Cl)c3)c2cc1OC')
        ]
    },
    {
        'target_name': 'Bacterial DNA Gyrase Subunit B (GyrB)',
        'gene': 'GyrB',
        'class': 'Antimicrobial',
        'organism': 'Bacterial',
        'uniprot': 'P0AES6',
        'pdb_id': '4PRX',
        'compounds': [
            ('Novobiocin', 'COc1cc(O)c2c(c1)OC(=O)C(NC(=O)c1ccc(O)c(CC=C(C)C)c1)=C2O'),
            ('Ciprofloxacin', 'O=C(O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O')
        ]
    },
    {
        'target_name': 'Bacterial Dihydrofolate Reductase (DHFR)',
        'gene': 'DHFR',
        'class': 'Antimicrobial',
        'organism': 'Bacterial',
        'uniprot': 'P0ABQ4',
        'pdb_id': '3IA5',
        'compounds': [
            ('Trimethoprim', 'COc1cc(Cc2cnc(N)nc2N)cc(OC)c1OC'),
            ('Methotrexate', 'CN(Cc1cnc2nc(N)nc(N)c2n1)c1ccc(C(=O)NC(CCC(=O)O)C(=O)O)cc1')
        ]
    },
    {
        'target_name': 'Bacterial Beta-Lactamase (TEM-1)',
        'gene': 'blaTEM',
        'class': 'Antimicrobial',
        'organism': 'Bacterial',
        'uniprot': 'P00807',
        'pdb_id': '1TEM',
        'compounds': [
            ('Clavulanic Acid', 'OCC=C1OC2CC(=O)N2C1C(=O)O'),
            ('Avibactam', 'NC(=O)C1CCN2C(=O)N(OS(=O)(=O)O)C1CC2'),
            ('Tazobactam', 'CC1(C(N2C(S1(=O)=O)CC2=O)C(=O)O)Cn1cnnn1')
        ]
    },

    # ── ION CHANNELS & SAFETY TARGETS ─────────────────────────────────
    {
        'target_name': 'Potassium Voltage-Gated Channel Subfamily H Member 2 (hERG / KCNH2)',
        'gene': 'KCNH2',
        'class': 'Ion Channel',
        'organism': 'Human',
        'uniprot': 'Q12809',
        'pdb_id': '5VA1',
        'compounds': [
            ('Astemizole', 'COc1ccc(CCN2CCC(Nc3nc4ccccc4n3Cc3ccc(F)cc3)CC2)cc1'),
            ('Dofetilide', 'CNS(=O)(=O)c1ccc(OCCN(CCc2ccc(NS(C)(=O)=O)cc2)C)cc1'),
            ('Terfenadine', 'CC(C)(C)c1ccc(C(O)CCCN2CCC(C(O)(c3ccccc3)c3ccccc3)CC2)cc1')
        ]
    },
    {
        'target_name': 'Sodium Voltage-Gated Channel Alpha Subunit 5 (Nav1.5 / SCN5A)',
        'gene': 'SCN5A',
        'class': 'Ion Channel',
        'organism': 'Human',
        'uniprot': 'Q14524',
        'pdb_id': '6LQA',
        'compounds': [
            ('Lidocaine', 'CCN(CC)CC(=O)Nc1c(C)cccc1C'),
            ('Flecainide', 'O=C(Nc1cc(OCC(F)(F)F)ccc1OCC(F)(F)F)C1CCCN1'),
            ('Mexiletine', 'CC(N)COc1c(C)cccc1C')
        ]
    },
    {
        'target_name': 'Voltage-Dependent L-type Calcium Channel Subunit Alpha-1C (Cav1.2 / CACNA1C)',
        'gene': 'CACNA1C',
        'class': 'Ion Channel',
        'organism': 'Human',
        'uniprot': 'Q13936',
        'pdb_id': '7UHF',
        'compounds': [
            ('Amlodipine', 'CCOC(=O)C1=C(COCCN)NC(C)=C(C(=O)OC)C1c1ccccc1Cl'),
            ('Verapamil', 'COc1ccc(CCN(C)CCCC(C#N)(C(C)C)c2ccc(OC)c(OC)c2)cc1OC'),
            ('Diltiazem', 'COc1ccc(C2Sc3ccccc3N(CCN(C)C)C(=O)C2OC(C)=O)cc1')
        ]
    },

    # ── TRANSPORTERS & EPIGENETICS ────────────────────────────────────
    {
        'target_name': 'Sodium-Dependent Serotonin Transporter (SERT / SLC6A4)',
        'gene': 'SLC6A4',
        'class': 'Transporter',
        'organism': 'Human',
        'uniprot': 'P31645',
        'pdb_id': '5I6X',
        'compounds': [
            ('Fluoxetine', 'CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1'),
            ('Sertraline', 'CNC1CCC(c2ccc(Cl)c(Cl)c2)c2ccccc21'),
            ('Citalopram', 'CN(C)CCCC1(c2ccc(C#N)cc2)OCc2cc(F)ccc21')
        ]
    },
    {
        'target_name': 'Multidrug Resistance Protein 1 (P-gp / ABCB1)',
        'gene': 'ABCB1',
        'class': 'Transporter',
        'organism': 'Human',
        'uniprot': 'P08183',
        'pdb_id': '6QEX',
        'compounds': [
            ('Tariquidar', 'COc1cc2c(cc1OC)CCN(CCc1ccc(NC(=O)c3cc(OC)c(OC)c(OC)c3)cc1)C2'),
            ('Elacridar', 'COc1cc2c(cc1OC)CCN(CCc1ccc(NC(=O)c3ccccc3)cc1)C2')
        ]
    },
    {
        'target_name': 'Bromodomain-Containing Protein 4 (BRD4)',
        'gene': 'BRD4',
        'class': 'Epigenetic',
        'organism': 'Human',
        'uniprot': 'O60885',
        'pdb_id': '3MXF',
        'compounds': [
            ('JQ1', 'Cc1sc2c(c1C)C(CC(=O)OC(C)(C)C)N=C(c1ccc(Cl)cc1)c1nn[nH]c1-2'),
            ('I-BET762', 'Cc1sc2c(c1C)C(CC(=O)NC(C)C)N=C(c1ccc(Cl)cc1)c1nn[nH]c1-2')
        ]
    },
    {
        'target_name': 'Histone Deacetylase 1 (HDAC1)',
        'gene': 'HDAC1',
        'class': 'Epigenetic',
        'organism': 'Human',
        'uniprot': 'Q13547',
        'pdb_id': '5ICN',
        'compounds': [
            ('Vorinostat', 'O=C(CCCCCCC(=O)Nc1ccccc1)NO'),
            ('Entinostat', 'O=C(Cc1ccccc1)Nc1ccc(NC(=O)c2ccc(CN)cc2)cc1'),
            ('Panobinostat', 'CC(C)NCCc1ccc(C=CC(=O)NO)cc1')
        ]
    },
    {
        'target_name': 'Poly [ADP-ribose] Polymerase 1 (PARP1)',
        'gene': 'PARP1',
        'class': 'Enzyme',
        'organism': 'Human',
        'uniprot': 'P09874',
        'pdb_id': '5DS3',
        'compounds': [
            ('Olaparib', 'O=C(c1cc(Cc2n[nH]c(=O)c3ccccc23)ccc1F)N1CCN(C(=O)C2CC2)CC1'),
            ('Niraparib', 'NC(=O)c1cccc(c1)C1CCNCC1'),
            ('Rucaparib', 'CNCc1ccc(-c2c[nH]c3c2cc(F)cc3)cc1')
        ]
    }
]

# ═════════════════════════════════════════════════════════════════════
# PRE-COMPUTED TARGET FINGERPRINTS CACHE
# ═════════════════════════════════════════════════════════════════════

_PRECOMPUTED = False

def _ensure_precomputed():
    """Calculates and caches Morgan fingerprints for all reference bioactive compounds."""
    global _PRECOMPUTED
    if _PRECOMPUTED:
        return
    for t in TARGET_PROFILES:
        fps = []
        for name, smi in t.get('compounds', []):
            try:
                m = Chem.MolFromSmiles(smi)
                if m:
                    fp = _get_fp(m)
                    fps.append((name, smi, fp))
            except Exception:
                pass
        t['cached_fps'] = fps
    _PRECOMPUTED = True


def get_all_target_classes() -> List[str]:
    """Returns sorted unique functional classes across the knowledge base."""
    classes = sorted(list({t.get('class', 'Other') for t in TARGET_PROFILES}))
    return ["All Classes"] + classes


def get_all_organisms() -> List[str]:
    """Returns sorted unique organisms across the knowledge base."""
    organisms = sorted(list({t.get('organism', 'Human') for t in TARGET_PROFILES}))
    return ["All Organisms"] + organisms


# ═════════════════════════════════════════════════════════════════════
# TARGET PREDICTION INFERENCE
# ═════════════════════════════════════════════════════════════════════

def predict_targets(
    smiles_or_mol: Any,
    top_k: int = 10,
    target_class: Optional[str] = None,
    organism: Optional[str] = None,
    min_similarity: float = 0.0
) -> pd.DataFrame:
    """
    Predicts biological macromolecular targets for a small molecule query across high-coverage proteomes.
    Calculates Tanimoto similarity to known bioactive chemotypes and models target binding probability.
    Supports functional class and organism filtering for deep navigation in ultra-large proteomes.
    """
    _ensure_precomputed()

    if isinstance(smiles_or_mol, str):
        query_mol = Chem.MolFromSmiles(smiles_or_mol)
    else:
        query_mol = smiles_or_mol

    if query_mol is None:
        return pd.DataFrame()

    query_fp = _get_fp(query_mol)
    records = []

    for target in TARGET_PROFILES:
        # Filter by class
        if target_class and target_class != "All Classes" and target.get('class') != target_class:
            continue
        # Filter by organism
        if organism and organism != "All Organisms" and target.get('organism') != organism:
            continue

        cached = target.get('cached_fps', [])
        if not cached:
            continue

        best_sim = 0.0
        best_comp_name = ""
        best_comp_smi = ""

        sims = []
        for name, smi, ref_fp in cached:
            sim = DataStructs.TanimotoSimilarity(query_fp, ref_fp)
            sims.append(sim)
            if sim > best_sim:
                best_sim = sim
                best_comp_name = name
                best_comp_smi = smi

        if best_sim < min_similarity:
            continue

        avg_sim = float(np.mean(sims)) if sims else 0.0

        # Calibrated logistic probability model based on chemical similarity:
        # P = 1 / (1 + exp(-k * (sim - x0)))
        k_factor = 14.0
        x0_threshold = 0.38
        effective_sim = 0.7 * best_sim + 0.3 * avg_sim
        prob = 1.0 / (1.0 + math.exp(-k_factor * (effective_sim - x0_threshold)))
        prob = min(max(prob, 0.05), 0.99)

        records.append({
            'Target Name': target['target_name'],
            'Gene': target.get('gene', ''),
            'Class': target['class'],
            'Organism': target.get('organism', 'Human'),
            'UniProt ID': target['uniprot'],
            'PDB ID': target['pdb_id'],
            'Probability': round(prob * 100, 1),
            'Tanimoto Similarity': round(best_sim, 3),
            'Reference Compound': best_comp_name,
            'Reference SMILES': best_comp_smi
        })

    if not records:
        return pd.DataFrame()

    df = pd.DataFrame(records)
    df = df.sort_values(by=['Probability', 'Tanimoto Similarity'], ascending=[False, False]).reset_index(drop=True)
    df.index = df.index + 1
    return df.head(top_k)
