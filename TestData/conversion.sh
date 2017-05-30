#!/bin/bash

for graphtype in gen2Mesh gen3Mesh genFixedLength genSegments
do
    #for weighttype in UnweightedUniformStretch WeightedUniformStretch UnweightedExpStretch WeightedExpStretch
    #do
    temp=${graphtype}WeightedExpStretch1e5
    python conversion.py A/${temp}.mtx T/${temp}_tree.mtx TplusE/${temp}_treeplus.01.mtx TplusE/${temp}_treeplus.05.mtx TplusE/${temp}_treeplus.1.mtx TplusE/${temp}_treeplus.15.mtx Factors/${temp}_factor.01.mtx Factors/${temp}_factor.05.mtx Factors/${temp}_factor.1.mtx Factors/${temp}_factor.15.mtx
    #done
done

