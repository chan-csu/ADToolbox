from Bio_Struct import *
from collections import Counter
import random
import csv
class Bio_seq:

    def __init__(self,seq="ATCG",seq_type="DNA",label='No Label'):
        self.seq=seq.upper()
        self.label=label
        self.seq_type=seq_type
        self.is_valid=self.__validateSeq()
        assert self.is_valid, f"Provided data does not seem to be a correct {self.seq_type} sequence!"
    def __validateSeq(self):
        return set(Nucleotides).issuperset(self.seq)
    def show_seq_info(self):
        return f"[Label]:{self.label}\n[Sequence]:{self.seq}\n[Molecule Type]:{self.seq_type}\n[Length]:{len(self.seq)}"
    def get_seq_biotype(self):
        return self.seq_type
    def generate_rand_seq(self,length=10,seq_type="DNA"):
        seq=''.join([random.choice(Nucleotides) for nuc in range(length)])
        self.__init__(seq,seq_type,"Randomly Generated DNA")

    def countNucFrequency(self):
        return dict(Counter(self.seq))

    def transcription(self):
        return self.seq.replace("T", "U")

    def reverse_complement(self):
       mapping=str.maketrans('ATCG','TAGC')
       return self.seq.translate(mapping)[::-1]

    def GC_Content(self):
        return round(((self.seq.count('C') + self.seq.count('G')) / len(self.seq) * 100), 6)

    def GC_Content_subseq(self, k=20):
        res = []
        for i in range(0, len(self.seq) - k + 1, k):
            subseq = self.seq[i:i + k]
            res.append(round(((subseq.count('C') + subseq.count('G')) / len(subseq) * 100), 6))
        return res

    def translate_seq(self, init_pos=0):
        return [DNA_Codons[self.seq[pos:pos + 3]] for pos in range(init_pos, len(self.seq) - 2, 3)]

    def codon_usage(self, aminoacid):
        tmpList = []
        for i in range(0, len(self.seq) - 2, 3):
            if DNA_Codons[self.seq[i:i + 3]] == aminoacid:
                tmpList.append(self.seq[i:i + 3])
        freqDict = dict(Counter(tmpList))
        totalweight = sum(freqDict.values())
        for seq in freqDict:
            freqDict[seq] = round(freqDict[seq] / totalweight, 2)
        return freqDict

    def gen_reading_frames(self):
        frames = []
        frames.append(self.translate_seq(0))
        frames.append(self.translate_seq(1))
        frames.append(self.translate_seq(2))
        tmp_seq=Bio_seq(self.reverse_complement(), self.seq_type)
        frames.append(tmp_seq.translate_seq(0))
        frames.append(tmp_seq.translate_seq(1))
        frames.append(tmp_seq.translate_seq(2))
        return frames

    def protein_from_rf(self,aa_seq):
        current_prot = []
        proteins = []
        for aa in aa_seq:
            if aa == "_":
                if current_prot:
                    for p in current_prot:
                        proteins.append(p)
                    current_prot = []
            else:
                if aa == "M":
                    current_prot.append("")
                for i in range(len(current_prot)):
                    current_prot[i] += aa
        return proteins

    def all_proteins_from_ORFs(self, startReadPos=0, endReadPos=0, ordered=False):

        if endReadPos > startReadPos:
            tmp_seq = Bio_seq(self.seq[startReadPos: endReadPos], self.seq_type)
            rfs = tmp_seq.gen_reading_frames()
        else:
            rfs = self.gen_reading_frames()
        res = []
        for rf in rfs:
            prots = self.protein_from_rf(rf)
            for p in prots:
                res.append(p)
        if ordered:
            return sorted(res, key=len, reverse=True)
        return res

    def RNA_To_Protein(self):
        Protein_seq=''
        if self.seq_type!="RNA":
            print("This function is only for RNA!")
        else:
            Sweep=0
            while Sweep+3<len(self.seq):
                Codon=self.seq[Sweep:Sweep+3]
                Protein_seq=Protein_seq+RNA_Codons[Codon]
                Sweep+=3

        return Protein_seq
    def Motif(self,Substring):
        K=len(Substring)
        L=len(self.seq)
        i=0
        Locations=''
        while (i+K)<L:
            if Substring==self.seq[i:i+K]:
                Locations+=' '+str(i+1)
            i+=1

        return Locations


class Fasta:

    def __init__(self,file):
        self.Fasta_File=file

    def Fasta_To_Dict(self):

        Dictionary={}

        for lines in self.Fasta_File:

            if '>' in lines:
                label = lines[1:-1]
                Dictionary[label]=""

            else:
                if '\n' in lines:
                    Dictionary[label] += lines[:-1]
                else:
                    Dictionary[label] += lines


        return Dictionary

    def overlap_graph(self,k):
        graph=[]
        Dict=self.Fasta_To_Dict()
        keylist =list(Dict.keys())
        for keys in keylist:
            for KEYS in keylist:
                if keys != KEYS:
                    if Dict[keys][-k:]==Dict[KEYS][0:k]:
                        graph.append((keys,KEYS))

        return graph



    def Shared_Motif(self):
        Dict=self.Fasta_To_Dict()
        Min_Length=min([len(Dict[key]) for key in Dict])
        Key_names=list(Dict.keys())
        Number_of_seqs=len(Key_names)
        Max_Substring='A'
        Sliding_Window=1
        while Sliding_Window<Min_Length:
            for i in range(len(Dict[Key_names[0]])-Sliding_Window):
                substring=Dict[Key_names[0]][i:i+Sliding_Window]
                j=0
                for keys in Key_names:
                    if Dict[keys].__contains__(substring):
                        j+=1
                if j==Number_of_seqs:
                    Max_Substring=substring
            Sliding_Window+=1



        return Max_Substring
