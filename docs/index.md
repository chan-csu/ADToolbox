# Welcome to ADToolbox

Parsa Ghadermazi 
parsa96@colostate.edu

AD Toolbox is developed in Chan Lab at Colorado State University. The main goal of this toolbox is to provide the tools that are useful for modeling and optimization of anaerobic digestion process. The following diagram shows the ADToolbox workflow:

```mermaid

flowchart LR
    mod4 & mod1 & pdb --> mod2
    mod2 & rdb --> mod3 
    mod3 --> ADM
    subgraph  MM[Metagenomics]
    style MM text-align:center
    mod1( Amplicon to Genome )
    mod2( Align Genomes )
    mod3( ADM Mapping )
    mod4( Input Genomes )
    mod0( 16s from QIIME II) 
    mod0 --> mod1

    end
    subgraph ADM[ ADM ]
    ModADM[Modified-ADM]
    OrADM[Original-ADM]
    ADMDash[ ADM Dashboard]
    CSVrep[CSV Report]
    Esch[Escher Map]
    ModADM & OrADM --> ADMDash & CSVrep & Esch
    end
    subgraph DB
    pdb[Protein Database]
    rdb[Reaction Database]
    MP[Model Parameters]
    BP[Base Parameters]
    IC[Initial Conditions]
    InC[Inlet Conditions]
    rxn[Reactions]
    sp[Species]
    end 
    subgraph Optimization
    PF[Parameter Fitting]
    PO[Process Optimization]

    end 
    subgraph Expd[Experimental Data]
    end
    MP & BP & IC & InC & rxn & sp --> ADM
    Expd --> PF
    PF --> MP
    PO --> BP & IC & InC

```




