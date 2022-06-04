#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pandas as pd
import os
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import DataStructs



df = pd.DataFrame()

target = os.listdir("target_file/")
ref = os.listdir("reference_file/")




ref_suppl = Chem.SDMolSupplier(str("reference_file/"+ref[0]))
ref_m1 = ref_suppl[0]
smi_ref_suppl = Chem.MolToSmiles(ref_m1)

for i in range(0,len(target)):
    suppl = Chem.SDMolSupplier(str("target_file/"+target[i]))
    m1 = suppl[0]
    smi_suppl = Chem.MolToSmiles(m1)
    ms = [Chem.MolFromSmiles(smi_ref_suppl), Chem.MolFromSmiles(smi_suppl)]
    fps = [Chem.RDKFingerprint(x) for x in ms]
    sim = DataStructs.FingerprintSimilarity(fps[0],fps[1])
    df[str(target[i])]=[sim]


df = df.transpose()

df.to_csv("tanimoto_matrix.csv")

