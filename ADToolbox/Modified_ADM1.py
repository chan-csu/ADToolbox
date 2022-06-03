import graphviz as gv

g = gv.Digraph('G', filename='Process_Graph.gv')
# g.node("A","Feed")
# g.node("B","AminoAcids")
# g.node("C","Sugars")
# g.node("D","Lipids")
# g.node("E","Soluble Inerts")
# g.node("F","Particulate Inerts")
# g.view()
g.attr(compound='true')
g.attr('node', shape='box')
Input=["Feed"]
First_nodes=["TDS_Degradation","TSS_Degradation"]
Second_nodes=["AminoAcids","Sugars","Lipids","Soluble Inerts","Particulate Inerts"]
Third_nodes=["Lactate","Ethanol","Acetate","Propionate"]
Fourth_nodes=["Butyrate","Valerate"]




g.edge("Feed","TDS_Degradation")
g.edge("Feed","TSS_Degradation")
g.edge("TDS_Degradation","AminoAcids")
g.edge("TSS_Degradation","AminoAcids")
g.edge("TDS_Degradation","Sugars")
g.edge("TSS_Degradation","Sugars")
g.edge("TDS_Degradation","Lipids")
g.edge("TSS_Degradation","Lipids")
g.edge("TDS_Degradation","Soluble Inerts")
g.edge("TSS_Degradation","Particulate Inerts")





g.edge("AminoAcids","Lactate")
g.edge("AminoAcids","Ethanol")
g.edge("AminoAcids","Acetate")
g.edge("Sugars","Lactate")
g.edge("Sugars","Ethanol")
g.edge("Sugars","Acetate")
g.edge("Lipids","Lactate")
g.edge("Lipids","Ethanol")
g.edge("Lipids","Acetate")
g.edge("Acetate","Butyrate")
g.edge("Ethanol","Butyrate")
g.edge("Lactate","Butyrate")
g.edge("Sugars","Propionate")
g.edge("Ethanol","Valerate")
g.edge("Lactate","Valerate")
g.edge("Butyrate","Caproate")
g.edge("Lactate","Caproate")
g.edge("Ethanol","Caproate")
g.edge("Propionate","Valerate")


with g.subgraph(name='Input') as s:
    s.attr(rank = 'n1',color='lightgray',style='filled',label='Process')
    for n in Input: s.node(n)
    s.edges=([("Feed","TDS_Degradation")])
    

with g.subgraph(name='Feed') as s:
    s.attr(rank = 'n2')
    for n in First_nodes: s.node(n)

with g.subgraph(name='digestate') as s:
    s.attr(rank = 'n3',color='blue',label='Digestate')
    for n in Second_nodes: s.node(n)

with g.subgraph(name='SCFAs') as s:
    s.attr(rank = 'n4',color='blue',label='SCFAs')
    for n in Third_nodes: s.node(n)

with g.subgraph(name='Products') as s:
    s.attr(rank = 'n5',color='blue',label='MCFAs')
    for n in Fourth_nodes: s.node(n)

with g.subgraph(name='Products') as s:
    s.attr(rank = 'n4',color='blue',label='MCFAs')
    for n in Fourth_nodes: s.node(n)


g.save("Process_Graph.gv")
