# ----------- Importing modules necessary to call the class object:

import h5py
import numpy as np
import mappy as mp

print("Loading reference genome...\n")
honeybee_genome = mp.Aligner("/Users/kristinaulicna/Documents/Rotation_1/GCF_000002195.4_Amel_4.5_genomic.fa.gz",
                           preset="map-ont")

class Alignment:
    def __init__(self, fast5_file, txt_file_directory):
        """ Insert input file in fast5 format (= HDF5 file) name into "" with a full directory (from /Users/...).
            Insert the directory (from /Users/...) where you would like the output file in txt format to be stored,
            - don't forget the '/' at the end of the desired directory."""
        input_file = h5py.File(fast5_file, "r")
        readSeq = np.array(input_file["Analyses"]["Basecall_1D_000"]["BaseCalled_template"]["Fastq"])
        readSeq = str(readSeq).split("\\n")
        self.readSeq = readSeq[1]
        self.readID = readSeq[0].split("@")[1].split("_")[0]
        self.events = np.array(input_file["Analyses"]["Basecall_1D_000"]["BaseCalled_template"]["Events"])
        fast5_name = str(input_file).split("_")
        input_file.close()
        self.fast5_name = fast5_name[-5] + "_" + fast5_name[-4] + "_" + fast5_name[-3] + "_" + fast5_name[-2]
        self.output_file = open(txt_file_directory + "/" + self.readID + ".txt", "w")

    def AlignQueToRef(self):
        """ Creates an alignment of query sequence to the reference contig.
            Keeps track of multiple signals per single contig positions.
            Skips insertions and deletions in the aligned sequences. """

        hit_counter = 0
        for self.hit in honeybee_genome.map(self.readSeq):
            hit_counter += 1
        if hit_counter == 0:
            print ("Warning:", "No hits identified")
            #raise Exception("Warning", "No hits identified")
        if hit_counter > 1:
            print ("Warning:", "Multiple hits identified")
            #raise Exception("Warning", "Multiple hits identified")
        query_seq = self.readSeq[self.hit.q_st:self.hit.q_en]
        self.reference = honeybee_genome.seq(self.hit.ctg, self.hit.r_st, self.hit.r_en)
        self.ref_plus = self.reference
        if self.hit.strand == -1:
            self.reference = self.reference[::-1]
            self.reference = str(self.reference) \
                .replace('T', '%temp1%').replace('A', 'T').replace('%temp1%', 'A') \
                .replace('C', '%temp2%').replace('G', 'C').replace('%temp2%', 'G')

        # Header: ["RefP"(contig), "QueP", "EveP", "Signal", "Shift", "RefB", "QueB", "KmerB", "RefKmer", "QueKmer":]
        self.my_list_1 = []
        position_ref = 0
        position_que = 0
        position_events = 1
        itrevents = iter(self.events)
        while position_events < self.hit.q_st:
            line = next(itrevents)
            if line[5] != 0:
                position_events += line[5]
        for parameter in self.hit.cigar:
            if parameter[1] == 0:
                for i in range(parameter[0]):
                    refBase = self.reference[position_ref]
                    queBase = query_seq[position_que]
                    if refBase == queBase:
                        while position_events <= (position_que + self.hit.q_st):
                            if position_events == (position_que + self.hit.q_st):
                                que_kmer = str(line[4])
                                if position_ref == 0:
                                    ref_kmer = '__' + self.reference[position_ref : position_ref + 3]
                                elif position_ref == 1:
                                    ref_kmer = '_' + self.reference[position_ref : position_ref + 4]
                                else:
                                    ref_kmer = self.reference[position_ref - 2 : position_ref + 3]
                                if len(ref_kmer) < 5:
                                    ref_kmer = ref_kmer + ("_" * (5 - len(ref_kmer)))
                                if refBase != que_kmer[4] or queBase != que_kmer[4] or refBase != queBase:
                                    print("MisMatch:", "RefBase, QueBase and K-merBase are not matching",
                                          position_events)
                                    #raise Exception("MisMatch", "RefBase, QueBase and K-merBase are not matching", position_events)
                                if self.hit.strand == -1:
                                    self.my_list_1.append([str(self.hit.r_en - position_ref), str(self.hit.q_st + position_que),
                                                    str(position_events), str(line[0]), str(line[5]),
                                                    str(refBase), str(queBase), que_kmer[4], ref_kmer, que_kmer[2:7]])
                                if self.hit.strand == 1:
                                    self.my_list_1.append([str(self.hit.r_st + position_ref), str(self.hit.q_st + position_que),
                                                    str(position_events), str(line[0]), str(line[5]),
                                                    str(refBase), str(queBase), que_kmer[4], ref_kmer, que_kmer[2:7]])
                            try:
                                line = next(itrevents)
                            except StopIteration:
                                break
                            if line[5] != 0:
                                position_events += line[5]
                    position_ref += 1
                    position_que += 1
            elif parameter[1] == 2:
                position_ref += parameter[0]
            elif parameter[1] == 1:
                position_que += parameter[0]
        return self.my_list_1

    def ExtractCpGsignal(self):
        """ Creates an extraction of signal data for 5 positions tracking the CpG dinucleotides movement.
            Produces averaged signals per respective k-mer, subtracted from the Nanopolish model k-mer signal. """

        # Map CpG positions as the dinucleotide moves along the strand.
        contig_lists = [[] for i in range(4)]
        try:
            for line in self.my_list_1:
                if line[8][0:2] == "CG":
                    contig_lists[0].append(int(line[0]))
                if line[8][1:3] == "CG":
                    contig_lists[1].append(int(line[0]))
                if line[8][2:4] == "CG":
                    contig_lists[2].append(int(line[0]))
                if line[8][3:5] == "CG":
                    contig_lists[3].append(int(line[0]))
            for small_list in contig_lists:
                list(set(small_list))
                sorted(small_list)
        except AttributeError:
            raise Exception("Warning", "Run AlignQueToRef() function first")

        # Distinguish between positive & negative strand processing:
        schematic_list = ["'_ _ _ _ C'", "'_ _ _ C G'", "'_ _ C G _'", "'_ C G _ _'", "'C G _ _ _'"]
        strand = ''
        if self.hit.strand == 1:
            strand = '+'
            normalizer = 0
            contig_lists = [x - 4 for x in contig_lists[0]] + [x - 3 for x in contig_lists[1]] + \
                           [x - 2 for x in contig_lists[2]] + [x - 1 for x in contig_lists[3]]
            contig_lists = list(set(contig_lists))
            contig_lists = sorted(contig_lists)
            contig_lists = [contig_lists,
                            [x + 1 for x in contig_lists], [x + 2 for x in contig_lists],
                            [x + 3 for x in contig_lists], [x + 4 for x in contig_lists]]
            ctg_pos_lists = contig_lists
        elif self.hit.strand == -1:
            strand = '-'
            normalizer = 1
            contig_lists = [x + 1 for x in contig_lists[3]] + [x + 2 for x in contig_lists[2]] + \
                           [x + 3 for x in contig_lists[1]] + [x + 4 for x in contig_lists[0]]
            contig_lists = list(set(contig_lists))
            contig_lists = sorted(contig_lists, reverse=True)
            contig_lists = [contig_lists,
                            [x - 1 for x in contig_lists], [x - 2 for x in contig_lists],
                            [x - 3 for x in contig_lists], [x - 4 for x in contig_lists]]
            ctg_pos_lists = contig_lists[::-1]
            schematic_list = sorted(schematic_list, reverse=True)
        CpG_positions = []
        for schematic, small_list in zip(schematic_list, contig_lists):
            CpG_positions.append([schematic, "\t", len(small_list), "\t", small_list])

        # Create a dictionary with CpG 'ctg_pos' as keys and [signal list] as values:
        dictionary = {}
        for small_list in contig_lists:
            for number in small_list:
                for event in self.my_list_1:
                    if str(number) == event[0]:
                        kmer = event[8]
                        if number not in dictionary:
                            dictionary[number] = [kmer, []]
                        if kmer != dictionary[number][0]:
                            print ("Warning", "K-mers not identical for same contig position", kmer, dictionary[number][0])
                            raise Exception("Warning", "K-mers not identical for same contig position", kmer, dictionary[number][0])
                        dictionary[number][1].append(float(event[3]))
                if number not in dictionary:
                    dictionary[number] = [None, None]

        # Calculate the signal mean & subtract from respective Nanopolish model k-mer signal:
        from Nanopolish import ModelKmerValues
        model_dict = ModelKmerValues()
        for key, value in dictionary.items():
            if type(dictionary[key][0]) == str:
                dictionary[key][1] = np.mean(value[1])
            if value[0] in model_dict.keys():
                dictionary[key][1] = float(value[1]) - float(model_dict[value[0]])

        # Header: ["Ctg", "CtgPos", "Strand", "ReadID", "Fast5file", "RepRefKmer[3rd CpG]", "Values:" "____C", "___CG", "__CG_", "_CG__", "CG___"]
        self.my_list_2 = []
        for CpG_middle in ctg_pos_lists[2]:
            rep_kmer = self.ref_plus[((CpG_middle - self.hit.r_st) - (2 + normalizer)) : ((CpG_middle - self.hit.r_st) + (3 - normalizer))]
            self.my_list_2.append([self.hit.ctg, CpG_middle - normalizer, strand, self.readID, self.fast5_name, rep_kmer,
                                dictionary[CpG_middle - 2][1], dictionary[CpG_middle - 1][1], dictionary[CpG_middle][1],
                                dictionary[CpG_middle + 1][1], dictionary[CpG_middle + 2][1]])
        return self.my_list_2

    def StoreFileOutput(self):
        """ Creates a new file (name = readID) and stores the output lists
            of both functions in it (my_list_1 and my_list_2), line by line."""

        try:
            string = ''
            for single_line in self.my_list_2:
                for parameter in single_line:
                    string += (str(parameter) + "\t")
                string = string[:-1]
                string += "\n"
            self.output_file.write(string)
            self.output_file.close()

        except AttributeError:
            raise Exception("Warning", "Run AlignQueToRef() & ExtractCpGsignal() functions first")