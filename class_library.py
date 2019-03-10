import collections
import datetime
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
import tempfile
from io import BytesIO
import base64
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

# translate AA to DNA-seq
# remove cut-sites
# back-translate DNA-seq to AA = assert = equal
# insert intron
# yield result

class Tables():

    def __init__( self, aa_seq = '' ):
        self._build_aa_lookup()
        self._get_genetic_code()
        self._build_codon_groups()
        self.allowed_characters_nt = set( 'ATGC' )
        assert sorted(self.allowed_characters_nt ) == sorted(IUPAC.unambiguous_dna.letters)
        #self.import_kazusa_codon_table()
        #self.convert_codon_counts_2_freq()
        self._build_re_lookup()

        self.aa_seq = aa_seq
        self.cDNA_seq = ''
        self.cDNA_seq_cleaned = ''
        self.cDNA_seq_plus_i = ''

        return

    def _build_aa_lookup( self ):
        self.AA_lookup_dict = {
            'Ala' : { '1-letter' : 'A', 'name' : 'Alanine' },
            'Arg' : { '1-letter' : 'R', 'name' : 'Arginine' },
            'Asn' : { '1-letter' : 'N', 'name' : 'Asparagine' },
            'Asp' : { '1-letter' : 'D', 'name' : 'Aspartic acid' },
            'Cys' : { '1-letter' : 'C', 'name' : 'Cysteine' },
            'Glu' : { '1-letter' : 'E', 'name' : 'Glutamic acid' },
            'Gln' : { '1-letter' : 'Q', 'name' : 'Glutamine' },
            'Gly' : { '1-letter' : 'G', 'name' : 'Glycine' },
            'His' : { '1-letter' : 'H', 'name' : 'Histidine' },
            'Ile' : { '1-letter' : 'I', 'name' : 'Isoleucine' },
            'Leu' : { '1-letter' : 'L', 'name' : 'Leucine' },
            'Lys' : { '1-letter' : 'K', 'name' : 'Lysine' },
            'Met' : { '1-letter' : 'M', 'name' : 'Methionine' },
            'Phe' : { '1-letter' : 'F', 'name' : 'Phenylalanine' },
            'Pro' : { '1-letter' : 'P', 'name' : 'Proline' },
            'Ser' : { '1-letter' : 'S', 'name' : 'Serine' },
            'Thr' : { '1-letter' : 'T', 'name' : 'Threonine' },
            'Trp' : { '1-letter' : 'W', 'name' : 'Tryptophan' },
            'Tyr' : { '1-letter' : 'Y', 'name' : 'Tyrosine' },
            'Val' : { '1-letter' : 'V', 'name' : 'Valine' },
            'Stop' : { '1-letter' : '*', 'name' : 'Stop' },
            }
        return

    def _get_genetic_code( self ):
        self.codon2aa = {
            'TTT':'Phe',
            'TCT':'Ser',
            'TAT':'Tyr',
            'TGT':'Cys',
            'TTC':'Phe',
            'TCC':'Ser',
            'TAC':'Tyr',
            'TGC':'Cys',
            'TTA':'Leu',
            'TCA':'Ser',
            'TAA':'Stop',
            'TGA':'Stop',
            'TTG':'Leu',
            'TCG':'Ser',
            'TAG':'Stop',
            'TGG':'Trp',
            'CTT':'Leu',
            'CCT':'Pro',
            'CAT':'His',
            'CGT':'Arg',
            'CTC':'Leu',
            'CCC':'Pro',
            'CAC':'His',
            'CGC':'Arg',
            'CTA':'Leu',
            'CCA':'Pro',
            'CAA':'Gln',
            'CGA':'Arg',
            'CTG':'Leu',
            'CCG':'Pro',
            'CAG':'Gln',
            'CGG':'Arg',
            'ATT':'Ile',
            'ACT':'Thr',
            'AAT':'Asn',
            'AGT':'Ser',
            'ATC':'Ile',
            'ACC':'Thr',
            'AAC':'Asn',
            'AGC':'Ser',
            'ATA':'Ile',
            'ACA':'Thr',
            'AAA':'Lys',
            'AGA':'Arg',
            'ATG':'Met',
            'ACG':'Thr',
            'AAG':'Lys',
            'AGG':'Arg',
            'GTT':'Val',
            'GCT':'Ala',
            'GAT':'Asp',
            'GGT':'Gly',
            'GTC':'Val',
            'GCC':'Ala',
            'GAC':'Asp',
            'GGC':'Gly',
            'GTA':'Val',
            'GCA':'Ala',
            'GAA':'Glu',
            'GGA':'Gly',
            'GTG':'Val',
            'GCG':'Ala',
            'GAG':'Glu',
            'GGG':'Gly'
            }
        assert { k:self.AA_lookup_dict[v]['1-letter'] for k, v in self.codon2aa.items() if k not in ['TGA', 'TAA', 'TAG'] } == CodonTable.unambiguous_dna_by_name["Standard"].forward_table
        return

    def _build_codon_groups( self ):
        self.aa2codons = {}
        for aa in sorted( set( self.codon2aa.values() ) ):
            self.aa2codons[ aa ] = []
            for codon, _AA in sorted( self.codon2aa.items() ):
                if aa == _AA:
                    self.aa2codons[ aa ].append( codon )
        # additionally one-letter-lookup-dict
        self.codon2aa_oneletter = { codon : self.AA_lookup_dict[aa]['1-letter'] for codon, aa in self.codon2aa.items() }
        return

    def convert_codon_counts_2_freq( self ):
        self.aa2codon2freq = collections.defaultdict( dict )
        for aa, codons in self.aa2codons.items():
            sum_of_counts = 0.0
            for codon in codons:
                sum_of_counts += self.codon2count[ codon ]
            for codon in codons:
                count = self.codon2count[ codon ]
                freq = float( count ) / sum_of_counts
                self.aa2codon2freq[ aa ][ codon ] = round( freq , 2 )
        self.aa2codon_freq_list = {}
        for aa in self.aa2codon2freq:
            l = sorted( [ ( freq, codon ) for codon, freq in self.aa2codon2freq[ aa ].items() ] )
            l.reverse()
            self.aa2codon_freq_list[ aa ] = l
            self.aa2codon_freq_list[ aa ] = [ ( codon, freq ) for freq, codon in l ]
        # build dict: one-letter-aa to most-freq-codon:
        self.aa_oneletter_2_mostfreq_codon_dict = { self.AA_lookup_dict[aa]['1-letter'] : codon_freq_list[0][0] for aa, codon_freq_list in self.aa2codon_freq_list.items() } # { 'F', TTC }
        return
    
    def _kazusa_codon_table(self):
        s = '''UUU  5.0(  2110)  UCU  4.7(  1992)  UAU  2.6(  1085)  UGU  1.4(   601)
        UUC 27.1( 11411)  UCC 16.1(  6782)  UAC 22.8(  9579)  UGC 13.1(  5498)
        UUA  0.6(   247)  UCA  3.2(  1348)  UAA  1.0(   441)  UGA  0.5(   227)
        UUG  4.0(  1673)  UCG 16.1(  6763)  UAG  0.4(   183)  UGG 13.2(  5559)

        CUU  4.4(  1869)  CCU  8.1(  3416)  CAU  2.2(   919)  CGU  4.9(  2071)
        CUC 13.0(  5480)  CCC 29.5( 12409)  CAC 17.2(  7252)  CGC 34.9( 14676)
        CUA  2.6(  1086)  CCA  5.1(  2124)  CAA  4.2(  1780)  CGA  2.0(   841)
        CUG 65.2( 27420)  CCG 20.7(  8684)  CAG 36.3( 15283)  CGG 11.2(  4711)

        AUU  8.0(  3360)  ACU  5.2(  2171)  AAU  2.8(  1157)  AGU  2.6(  1089)
        AUC 26.6( 11200)  ACC 27.7( 11663)  AAC 28.5( 11977)  AGC 22.8(  9590)
        AUA  1.1(   443)  ACA  4.1(  1713)  AAA  2.4(  1028)  AGA  0.7(   287)
        AUG 25.7( 10796)  ACG 15.9(  6684)  AAG 43.3( 18212)  AGG  2.7(  1150)

        GUU  5.1(  2158)  GCU 16.7(  7030)  GAU  6.7(  2805)  GGU  9.5(  3984)
        GUC 15.4(  6496)  GCC 54.6( 22960)  GAC 41.7( 17519)  GGC 62.0( 26064)
        GUA  2.0(   857)  GCA 10.6(  4467)  GAA  2.8(  1172)  GGA  5.0(  2084)
        GUG 46.5( 19558)  GCG 44.4( 18688)  GAG 53.5( 22486)  GGG  9.7(  4087)'''
        return s
    
    def _hivecut_codon_table(self):
        s = '''TTT  5.41 (35482)   TCT  4.54 (29793)   TAT  2.77 (18152)   TGT  2.62 (17167) 
        TTC 20.32 (133224)  TCC 14.63 (95935)   TAC 18.40 (120642)  TGC 13.63 (89359) 
        TTA  0.84 (5538)    TCA  4.95 (32455)   TAA  0.45 (2957)    TGA  0.79 (5200)  
        TTG  6.33 (41523)   TCG 15.05 (98674)   TAG  0.53 (3475)    TGG 13.61 (89197) 

        CTT  5.36 (35150)   CCT  7.03 (46104)   CAT  3.46 (22676)   CGT  5.38 (35290) 
        CTC 14.13 (92666)   CCC 23.73 (155576)  CAC 18.38 (120500)  CGC 31.16 (204288)
        CTA  3.84 (25195)   CCA  7.24 (47485)   CAA  5.24 (34385)   CGA  3.35 (21931) 
        CTG 65.52 (429567)  CCG 26.94 (176646)  CAG 37.57 (246307)  CGG 19.05 (124861)

        ATT  5.82 (38140)   ACT  4.90 (32106)   AAT  3.13 (20520)   AGT  4.59 (30114) 
        ATC 18.05 (118325)  ACC 20.80 (136358)  AAC 20.47 (134183)  AGC 27.25 (178637)
        ATA  1.93 (12644)   ACA  6.57 (43096)   AAA  3.08 (20203)   AGA  1.16 (7596)  
        ATG 21.30 (139614)  ACG 17.80 (116694)  AAG 28.40 (186214)  AGG  4.97 (32581) 

        GTT  4.95 (32478)   GCT 17.60 (115362)  GAT  8.02 (52556)   GGT 10.50 (68847) 
        GTC 13.68 (89717)   GCC 51.11 (335087)  GAC 38.09 (249718)  GGC 61.79 (405085)
        GTA  2.94 (19269)   GCA 18.41 (120724)  GAA  4.49 (29453)   GGA  7.02 (46029) 
        GTG 46.65 (305828)  GCG 58.72 (384988)  GAG 46.98 (307968)  GGG 16.53 (108348)'''
        return s

    def import_codon_table( self, codon_table_string ):
        self.codon2count = {}
        codon_table_string = codon_table_string.replace( 'U', 'T' )
        for line in codon_table_string.split( '\n' ):
            line = line.strip()
            if line == [ '' ]:
                continue
            for codon_string in line.split( ')' )[ : -1 ]:
                codon_string = codon_string.strip()
                codon_freq, number = codon_string.split( '(' )
                number = int( number.strip() )
                codon, frequency_per_thousand = codon_freq.split()
                frequency_per_thousand = float( frequency_per_thousand.strip() )
                codon = codon.strip().upper()
                if self.check_codon( codon ):
                    self.codon2count[ codon.strip() ] = number
        return

    def check_codon( self, codon ):
        if len( codon ) != 3:
            return False
        if not set( codon ).issubset( self.allowed_characters_nt ):
            return False
        return True

    def _build_re_lookup( self ):
        self.re_lookup_dict = {
            'AatII' :  'GACGT^C',
            'AgeI'  :  'A^CCGGT',
            'AvrII' :  'C^CTAGG',
            'BamHI' :  'G^GATCC',
            'BbsI'  :  'GAAGAC',
            'BglII' :  'A^GATCT',
            'BsaI'  :  'GGTCTC',
            'EcoRI' :  'G^AATTC',
            'EcoRV' :  'GAT^ATC',
            'HindIII': 'A^AGCTT',
            'KpnI'  :  'GGTAC^C',
            'MluI'  :  'A^CGCGT',
            'NdeI'  :  'CA^TATG',
            'ScaI'  :  'AGT^ACT',
            'SmaI'  :  'CCC^GGG',
            'SpeI'  :  'A^CTAGT',
            'XbaI'  :  'T^CTAGA',
            'XhoI'  :  'C^TCGAG',
            'ZraI'  :  'GAC^GTC'
            }
        self.re_lookup_dict = { k : v.replace('^', '') for k, v in self.re_lookup_dict.items() }

        return

    def parse_fasta( self, fasta_seq ):
        seq_dict = collections.OrderedDict()
        name = 'unnamed_input_seq'
        newline = '{0}{1}'.format( chr(13), chr(10) )
        if not newline in fasta_seq:
            newline = '{0}'.format( chr(10) )
        if not newline in fasta_seq:
            newline = '{0}'.format( chr(13) )
        for line in fasta_seq.split( newline ):
            if line.startswith('>'):
                name = line[1:].strip() #remove the arrow >
                #name = name.replace( '_', ' ' )
                seq_dict[ name ] = ''
            else:
                try:
                    seq_dict[ name ] += line.strip().upper()
                except KeyError:
                    seq_dict[ name ] = line.strip().upper()
        return seq_dict



class CodonOptimizer():

    def __init__( self, aa_seq ):
        self.aa_seq = aa_seq
        self.dna_seq = ''
        return

    def reverse_translate( self, aa2codon_dict ):
        for aa in self.aa_seq:
            self.dna_seq += aa2codon_dict[ aa.upper() ]
            # what if invalid characters?
        return

    def translate( self, dna_seq, codon2aa_dict ):
        codon = ''
        translated_aa_seq = ''
        for i, n in enumerate( dna_seq.upper() ):
            i += 1
            codon += n
            if i % 3 == 0:
                translated_aa_seq += codon2aa_dict[ codon ] # what if invalid codon?
                codon = ''
        return translated_aa_seq

    def reverse_complement( self, dna_seq ):
        lookup_dict = { 'A' : 'T', 'T' : 'A', 'G' : 'C', 'C' : 'G' }
        dna_seq_new = ''
        for character in reversed( dna_seq.upper() ):
            dna_seq_new += lookup_dict[ character ]
        return dna_seq_new



class CutSiteRemover():

    def __init__( self, dna_seq, cut_site_list, codon2aa, aa2codon_freq_list ):
        self.messages = []
        self.dna_seq = dna_seq.upper()
        self.dna_seq_new = dna_seq
        self.cut_site_list = [ cut_site.upper() for cut_site in cut_site_list ]
        # add reverse complement of cut sites to list
        CO_class = CodonOptimizer( aa_seq = '' )
        self.cut_site_list += [ CO_class.reverse_complement( nt_seq ) for nt_seq in self.cut_site_list ]
        self.cut_site_list = sorted( set( self.cut_site_list ) )
        self.codon2aa = codon2aa
        self.aa2codon_freq_list = aa2codon_freq_list
        if self.cut_site_list:
            self.max_cut_site_len = max([len(c) for c in self.cut_site_list])
        else:
            self.max_cut_site_len = None
        self.log_dict = {}
        return
    
    def substring_indexes( self, substring, string ):
        last_found = -1  # Begin at -1 so the next position to search from is 0
        while True:
            # Find next index of substring, by starting after its last known position
            last_found = string.find(substring, last_found + 1)
            if last_found == -1:
                break  # All occurrences have been found
            yield last_found

    def get_cut_site_indices( self, dna_seq = '' ):
        if dna_seq == '':
            dna_seq = self.dna_seq
        dna_seq = dna_seq.upper()
        cut_site2index_list = collections.defaultdict( list )
        for cut_site in self.cut_site_list:
            for last_found_index in self.substring_indexes( substring = cut_site, string = dna_seq ):
                cut_site2index_list[ cut_site ].append( last_found_index )
        return cut_site2index_list

    def determine_frame( self, index_cut_site, cut_site ):
        len_cut_site = len( cut_site )
        i = index_cut_site
        # the following can surely be transferred to 10er cutter etc. as well, maybe by
        # length_of_affected_dna_seq (12 or 9 or 6) - length_cut_site (8 or 6) - start_shift (0 or 1 or 2) = addition (8er: 12-8-2 = 2 or 9-8-1 = 0)
        if len_cut_site == 6:
            if i % 3 == 0: # in frame
                index_2_start_from = i
                addition = 0
            elif (i - 1) % 3 == 0: # frame starts at -1
                index_2_start_from = i - 1
                addition = 2
            elif (i - 2) % 3 == 0: # frame starts at -2
                index_2_start_from = i - 2
                addition = 1
        elif len_cut_site == 8:
            if i % 3 == 0: # 3 codons can be considered for replacement
                index_2_start_from = i
                addition = 1
            elif (i - 1) % 3 == 0: # 3 codons can be considered for replacement
                index_2_start_from = i - 1
                addition = 0
            elif (i - 2) % 3 == 0: # 4 codons can be considered for replacement
                index_2_start_from = i - 2
                addition = 2
        else:
            self.messages.append( '[ ERROR ] the length of the cut site {0} is not supported, i. e. has not 6 or 8 nt characters.'.format( cut_site ) )
        return index_2_start_from, addition

    def get_codon_list( self, dna_seq ):
        if len(dna_seq) % 3 != 0:
            raise Exception('[ ERROR ] dna_seq "{0}" is not a multiplier of 3! Thus codons cannot be determined...'.format( dna_seq ))
        codon_list = []
        codon = ''
        for j, nt in enumerate( dna_seq ):
            j += 1
            codon += nt
            if j % 3 == 0:
                codon_list.append( codon )
                codon = ''
        return codon_list

    def check_codon_replacement_removes_cutsite ( self, old_codon, new_codon, codon_list, cut_site ):
        dna_seq = ''.join( codon_list )
        if dna_seq.find( cut_site ) == -1:
            raise Exception('[ ERROR ] current dna seq "{0}" does not contain cut site "{1}".'.format( dna_seq, cut_site ))
        dna_seq_new = self.replace_codon(codon_list, old_codon, new_codon)
        if dna_seq_new.find( cut_site ) == -1: # cut site removed
            return True
        else: # cut site not removed
            return False

    def replace_codon(self, codon_list, old_codon, new_codon):
        codon_list_copy = list(codon_list)
        pos = codon_list_copy.index(old_codon)
        codon_list_copy[pos] = new_codon
        return ''.join(codon_list_copy)

    # update 27.02.2019
    def check_replacement_introduces_other_cut_site(self, codon_index, index_2_start_from, dna_seq):
        position_replaced_codon = index_2_start_from + codon_index * 3
        # the following only works for cut sites with 6 nt - not 8er cutter!
        frame_min = position_replaced_codon - (self.max_cut_site_len - 1) # either 5 or 7, for 6er or 8er cutter
        if frame_min < 0:
            frame_min = 0
        frame_max = position_replaced_codon + (self.max_cut_site_len - 1) + 3
        if frame_max > len(dna_seq):
            frame_max = len(dna_seq)
        dna_seq_original = self.dna_seq_new[frame_min : frame_max]
        dna_seq_cleaned = dna_seq[frame_min : frame_max]
        n_cut_sites_original, n_cut_sites_cleaned = 0, 0
        for cut_site in self.cut_site_list:
            for i in self.substring_indexes(substring=cut_site, string=dna_seq_original):
                if i != -1:
                    n_cut_sites_original += 1
        for cut_site in self.cut_site_list:
            for i in self.substring_indexes(substring=cut_site, string=dna_seq_cleaned):
                if i != -1:
                    n_cut_sites_cleaned += 1
        #print(index_2_start_from, codon_index, dna_seq_original, dna_seq_cleaned, n_cut_sites_original, n_cut_sites_cleaned)
        if n_cut_sites_cleaned < n_cut_sites_original:
            return True
        return False


    def main( self, iter_max = 1000 ):
        iter_counter = 0
        self.log_dict = {}
        while iter_counter < iter_max:
            iter_counter += 1
            cut_site_found = False
            for cut_site in self.cut_site_list:
                for i in self.substring_indexes( substring = cut_site, string = self.dna_seq_new ):
                    if cut_site not in self.log_dict:
                        self.log_dict[ cut_site ] = {}
                    self.log_dict[ cut_site ][ i ] = {}
                    self.log_dict[ cut_site ][ i ][ 'freq_altcodon_orgcodon_remove_list' ] = []
                    cut_site_found = True
                    # get frame
                    index_2_start_from, addition = self.determine_frame( index_cut_site = i, cut_site = cut_site )
                    self.log_dict[ cut_site ][ i ][ 'index_2_start_from' ] = index_2_start_from
                    self.log_dict[ cut_site ][ i ][ 'region_with_cut_site' ] = '{0}'.format( self.dna_seq_new[ index_2_start_from : i + len( cut_site ) + addition ] )
                    # get codons to replace
                    codon_list = self.get_codon_list( dna_seq = self.dna_seq_new[ index_2_start_from : i + len( cut_site ) + addition ] )
                    freq_codon_list = []
                    for j, codon in enumerate(codon_list):
                        codon_freq_list = self.aa2codon_freq_list[ self.codon2aa[ codon ] ]
                        if len( codon_freq_list ) == 1:
                            self.log_dict[ cut_site ][ i ][ 'freq_altcodon_orgcodon_remove_list' ].append( ( codon_freq_list[0][-1], codon_freq_list[0][0], codon_freq_list[0][0], False ) )
                            continue # codon cannot be exchanged = ATG or TGG, thus go on
                        old_codon_index = [ c for c, f in codon_freq_list ].index( codon )
                        old_codon, old_freq = codon_freq_list[ old_codon_index ]
                        if codon != old_codon:
                            raise Exception('[ ERROR ] current codon = "{0}", old codon = "{1}"'.format( codon, old_codon ))
                        self.log_dict[ cut_site ][ i ][ 'freq_altcodon_orgcodon_remove_list' ].append( ( old_freq, old_codon, old_codon, False ) )
                        for c, f in codon_freq_list[ old_codon_index + 1 : ]:
                            if not self.check_codon_replacement_removes_cutsite(old_codon=codon, new_codon=c, codon_list=codon_list, cut_site=cut_site):
                                self.log_dict[ cut_site ][ i ][ 'freq_altcodon_orgcodon_remove_list' ].append( ( f, c, old_codon, False ) )
                                continue # codon replacement would not remove the cut site, thus go on
                            # update 27.02.2019
                            dna_seq_cut_site_removed = self.replace_codon(codon_list, old_codon, new_codon=c)
                            if not self.check_replacement_introduces_other_cut_site(
                                    codon_index=j,
                                    index_2_start_from=index_2_start_from,
                                    dna_seq=self.dna_seq_new[ : index_2_start_from ] + dna_seq_cut_site_removed + self.dna_seq_new[ i + len( cut_site ) + addition : ]):
                                #self.log_dict[ cut_site ][ i ][ 'freq_altcodon_orgcodon_remove_list' ].append( ( f, c, old_codon, 'False2' ) )
                                self.log_dict[ cut_site ][ i ][ 'freq_altcodon_orgcodon_remove_list' ].append( ( f, c, old_codon, False ) )
                                continue
                            self.log_dict[ cut_site ][ i ][ 'freq_altcodon_orgcodon_remove_list' ].append( ( f, c, old_codon, True ) )
                            freq_codon_list.append( ( f, c, old_codon ) )
                    # determine codon to replace
                    if not freq_codon_list:
                        self.messages.append( '[ ERROR ] cut site {0}, found at {1}, frame starting at {2}, could not be removed! Please check manually.'.format( cut_site, i, index_2_start_from ) )
                        self.log_dict[ cut_site ][ i ][ 'new_codon' ] = 'ERROR - could not be removed'
                        self.log_dict[ cut_site ][ i ][ 'region_without_cut_site' ] = 'ERROR - could not be removed'
                    else:
                        max_new_freq, new_codon, old_codon = sorted( freq_codon_list)[ -1 ]
                        assert self.check_codon_replacement_removes_cutsite(old_codon=old_codon, new_codon=new_codon, codon_list=codon_list, cut_site=cut_site)
                        self.log_dict[ cut_site ][ i ][ 'new_codon' ] = new_codon
                        # insert new codon
                        dna_seq_cut_site_removed = self.replace_codon(codon_list, old_codon, new_codon)
                        assert dna_seq_cut_site_removed.find( cut_site ) == -1
                        # update dna_seq
                        self.dna_seq_new = self.dna_seq_new[ : index_2_start_from ] + dna_seq_cut_site_removed + self.dna_seq_new[ i + len( cut_site ) + addition : ]
                        self.log_dict[ cut_site ][ i ][ 'region_without_cut_site' ] = '{0}'.format( self.dna_seq_new[ index_2_start_from : i + len( cut_site ) + addition ] )
            if not cut_site_found:
                break
        if iter_counter == iter_max:
            self.messages.append( '[ ERROR ] cut sites could not be removed despite {0} iterations.'.format(iter_counter) )
        return self.dna_seq_new
    

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class IntronError(Error):
    """Exception raised for errors in the intron insertion process."""
    pass


        
class IntronInserter():

    def __init__( self, dna_seq, intron_seq, start = 100, intermediate = 450, end = 100, max = 500, insertion_seq = 'GG', start_position_list = None, end_intron_different_seq = False, CSR_class = None ):
        self.messages = []
        self.dna_seq = dna_seq.upper()
        self.insertion_seq = insertion_seq
        self.start = start
        self.intermediate = intermediate
        self.end = end
        self.max = max
        self.iter_max = 1000
        self.start_position_list = start_position_list
        if self.start_position_list == None:
            self.start_position_list = []
        self.position_list = []
        self.dna_seq_new = ''
        self.log_dict = {}
        self.intron_seq = intron_seq.lower()
        self.end_intron_different_seq = end_intron_different_seq.lower()
        if self.end_intron_different_seq == '':
            self.end_intron_different_seq = self.intron_seq
        self.CSR_class = CSR_class
        return


    def _determine_positions_from_manual_input(self):
        position_list = []
        try:
            for start_position in sorted( self.start_position_list ):
                possible_positions = self._get_insertion_positions( start_position = start_position, directions = [-1, 1], max_search_range = self.max )
                if not possible_positions:
                    self.messages.append( '[ ERROR ] Could not find the nucleotide pair "{2}" within -{0} nt around the position {1} during the process of inserting introns. No introns will be inserted. Please supply a refined position list or choose automatic insertion!'.format( self.max, start_position, self.insertion_seq ) )
                    raise IntronError
                for z in range(1000):
                    try:
                        position = possible_positions.pop()
                    except IndexError:
                        self.messages.append( '[ ERROR ] Inserting an intron was not possible. As the original position would have introduced a cut site, alternative positions were tried; however, not enough alternative positions were located in the sequence. Please supply a different position list or use the option "automatic determination". No introns will be inserted!' )
                        raise IntronError
                    cut_site2index_list = self._is_cut_site_introduced(position = position, intron_seq = self.intron_seq, CSR_class = self.CSR_class)
                    if cut_site2index_list:
                        self.messages.append( '[ INFO ] Inserting an intron at the original position {0} would introduce the cut site {1}, therefore an alternative position is tried ...'.format(position, cut_site2index_list) )
                    else:
                        break
                position_list.append( position )
            exon_length_list = [ position_list[0] + 1 ] + [ j - i for i, j in zip( position_list[:-1], position_list[1:] ) ] + [ len(self.dna_seq) - position_list[-1] - 1 ]
            if max( exon_length_list ) > self.max:
                self.messages.append( '[ WARNING ] The max exon length of {0} nt is exceeded in your supplied position list of {1}, which maps to the position list of {2} using the insertion seq "{3}", resulting into exon lengths of {4}. Better supply a refined position list to stay below {0} nt max exon length!'.format( self.max, self.start_position_list, position_list, self.insertion_seq, exon_length_list ) )
        except IntronError:
            position_list = []
        return position_list



    def determine_positions( self ):
        if self.start_position_list: 
            self.position_list = self._determine_positions_from_manual_input()
        else: # automatically_determine_intron_positions = True
            try:
                self.start_position_list = []
                self.log_dict[ 'seq_len' ] = len( self.dna_seq ) # alters the output in the LogClass
                self.log_dict[ 'seq_len_remaining' ] = None
                self.log_dict[ 'target_intermediate_exon_len' ] = self.intermediate
                self.log_dict[ 'initial_exon_determination' ] = None # [ min exon number, dist_1, max exon number, dist_2, chosen number of exons from this step ]
                self.log_dict[ 'final_exon_determination' ] = None # [ had the exon counter to be increased?, intermediate exon len, finally chosen exon number ]
                # insert end intron
                if len( self.dna_seq ) >= self.end: # 100
                    putative_position, max_search_range = len( self.dna_seq ) - self.end - 1, self.end - 2 # -1, because we are in positions = len(x) will raise IndexError
                    end_position_list = self._get_insertion_positions( start_position = putative_position, directions = [1], max_search_range = max_search_range )
                    if not end_position_list:
                        self.messages.append( '[ ERROR ] Unable to find the nucleotide pair "{2}" within +{0} nt upstream of the position {1}. Please increase the value for the parameter "end", which is currently set to {3}. No introns will be inserted!'.format( max_search_range, putative_position, self.insertion_seq, self.end ) )
                        raise IntronError
                    for z in range(1000): # instead of "while True:" to avoid being trapped in an endless loop
                        try:
                            end = end_position_list.pop()
                        except IndexError:
                            self.messages.append( '[ ERROR ] Inserting the "end-intron" was not possible. As the original position would have introduced a cut site, alternative positions were tried; however, not enough alternative positions were located in the sequence. Please increase the value for the parameter "end". No introns will be inserted!' )
                            raise IntronError
                        cut_site2index_list = self._is_cut_site_introduced(position = end, intron_seq = self.end_intron_different_seq, CSR_class = self.CSR_class)
                        if cut_site2index_list:
                            self.messages.append( '[ INFO ] Inserting the "end-intron" at the original "end-position" {0} would introduce the cut site {1}, therefore an alternative position is tried ...'.format(end, cut_site2index_list) )
                        else:
                            break
                else: # do not insert any intron if the seq is less than self.end, e.g. < 100 bp
                    raise IntronError
                # insert start intron
                if len( self.dna_seq ) >= self.start: # 100
                    putative_position, max_search_range = self.start - 1, self.start - 2 # -1, because GG searched and found at 100 and 101 would result in 101 nt start exon len
                    start_position_list = self._get_insertion_positions( start_position = putative_position, directions = [-1], max_search_range = max_search_range )
                    if not start_position_list:
                        self.messages.append( '[ ERROR ] Unable to find the nucleotide pair "{2}" within -{0} nt downstream of the position {1}. Please increase the value for the parameter "start", which is currently set to {3}. No introns will be inserted!'.format( max_search_range, putative_position, self.insertion_seq, self.start ) ) 
                        raise IntronError
                    for z in range(1000): # instead of "while True:" to avoid being trapped in an endless loop
                        try:
                            start = start_position_list.pop()
                        except IndexError:
                            self.messages.append( '[ ERROR ] Inserting the "start-intron" was not possible. As the original position would have introduced a cut site, alternative positions were tried; however, not enough alternative positions were located in the sequence. Please increase the value for the parameter "start". No introns will be inserted!' )
                            raise IntronError
                        cut_site2index_list = self._is_cut_site_introduced(position = start, intron_seq = self.intron_seq, CSR_class = self.CSR_class)
                        if cut_site2index_list:
                            self.messages.append( '[ INFO ] Inserting the "start-intron" at the original "start-position" {0} would introduce the cut site {1}, therefore an alternative position is tried ...'.format(start, cut_site2index_list) )
                        else:
                            break                    
                else: # do not insert the start intron if the seq is less than self.start, e.g. < 100 bp (self.end + self.start = 200, e.g. 199 bp => only end intron at pos 150 => 150 bp without intron, although start = 100 was requested)
                    self.start_position_list = [ len( self.dna_seq ) - self.end ]
                    self.position_list = [ end ]
                    raise IntronError
                # insert intermediate introns
                if len( self.dna_seq ) <= start + self.max + (len( self.dna_seq ) - end): # 700, if start=self.start and end=start.end
                    self.start_position_list = [ self.start, len( self.dna_seq ) - self.end ]
                    self.position_list = [ start, end ]
                    raise IntronError
                else: # thus > 700
                    start_position_list = []
                    remaining_len = float( end - start )
                    self.log_dict[ 'seq_len_remaining' ] = int( remaining_len )
                    min_exon_number = float( int( remaining_len / float( self.intermediate ) ) )
                    max_exon_number = min_exon_number + 1.0
                    dist_1 = abs( self.intermediate - remaining_len / min_exon_number )
                    dist_2 = abs( self.intermediate - remaining_len / max_exon_number )
                    if dist_1 < dist_2:
                        number_of_exons = min_exon_number
                    elif dist_1 > dist_2:
                        number_of_exons = max_exon_number
                    else: #equal
                        number_of_exons = min_exon_number
                    self.log_dict[ 'initial_exon_determination' ] = [ int(min_exon_number), int(dist_1), int(max_exon_number), int(dist_2), int(number_of_exons) ]
                    self.log_dict[ 'final_exon_determination' ] = ( True, None, None )
                    positions_satisfyingly_determined = False
                    for k in range( 1000 ):
                        if positions_satisfyingly_determined:
                            break
                        start_position_list, intermediate_exon_len = self._determine_intermediate_positions( number_of_exons = number_of_exons, start = start, remaining_len = remaining_len ) # empty, if number_of_exons == 1
                        if self.log_dict[ 'final_exon_determination' ][0] == True:
                            self.log_dict[ 'final_exon_determination' ] = ( True, intermediate_exon_len, int(number_of_exons) )
                        else:
                            self.log_dict[ 'final_exon_determination' ] = (False, intermediate_exon_len, int(number_of_exons))
                        max_search_range = intermediate_exon_len
                        possible_positions_list = []
                        for start_position in start_position_list:
                            position_list = self._get_insertion_positions( start_position = start_position, directions = [-1, 1], max_search_range = max_search_range )
                            if not position_list:
                                self.messages.append( '[ ERROR ] Could not find the nucleotide pair "{2}" within -{0} nt and +{0} nt around the position {1} during the process of inserting "intermediate-introns" ({3} exons). Please change the nucleotide pair or increase max exon length, which are currently set to "{2}" and {4} nt, respectively. No introns will be inserted!'.format( max_search_range, start_position, self.insertion_seq, number_of_exons, self.max ) )
                                raise IntronError
                            possible_positions_list.append( position_list )
                        # check if all exons remain below self.max (or are equal)
                        if len( possible_positions_list ) == 0: # possible_positions_list might be empty, if one exon!
                            exon_length_list = [ remaining_len ]
                            if max( exon_length_list ) <= self.max:
                                positions_satisfyingly_determined = True
                            else:
                                number_of_exons += 1
                                self.log_dict[ 'final_exon_determination' ] = ( False, intermediate_exon_len, number_of_exons )
                        else:
                            actual_position_list = [ l.pop() for l in possible_positions_list ] # l contains at least one element
                            for z in range(1000): # instead of "While True:" to avoid being trapped in an endless loop
                                exon_length_list = [ actual_position_list[0] - start ] + [ j - i for i, j in zip( actual_position_list[:-1], actual_position_list[1:] ) ] + [ end - actual_position_list[-1] ]
                                # if max distance between intermediate exons > self.max, distribute one more intermediate intron:
                                if max( exon_length_list ) <= self.max: # all exon lengths remain lower (or are equal) than max exon length
                                    position_would_introduce_cut_site = False
                                    for i, pos in enumerate(list(actual_position_list)): # check whether cut sites are introduced due to intron insertion
                                        cut_site2index_list = self._is_cut_site_introduced(position = pos, intron_seq = self.intron_seq, CSR_class = self.CSR_class)
                                        if cut_site2index_list: # if cut site found, try alternative position
                                            position_would_introduce_cut_site = True
                                            try:
                                                actual_position_list[i] = possible_positions_list[i].pop() # replace position with alternative position
                                            except IndexError:
                                                self.messages.append( '[ ERROR ] Inserting an "intermediate-intron" was not possible. As the original position would have introduced a cut site, alternative positions were tried; however, not enough alternative positions were located in the sequence. Please increase the value for the parameter "max exon length". No introns will be inserted!' )
                                                raise IntronError
                                            self.messages.append( '[ INFO ] Inserting an "intermediate-intron" at the original "intermediate-position" {0} would introduce the cut site {1}, therefore try alternative position ...'.format(pos, cut_site2index_list) )
                                    if not position_would_introduce_cut_site:
                                        positions_satisfyingly_determined = True
                                        break
                                else:
                                    number_of_exons += 1
                                    self.log_dict[ 'final_exon_determination' ] = ( False, intermediate_exon_len, number_of_exons )
                                    break
                    self.start_position_list = [ self.start ] + start_position_list + [ len( self.dna_seq ) - self.end ]
                    self.position_list = [ start ] + actual_position_list + [ end ]
            except IntronError:
                pass
        # sort, because intron insertion logic = self.insert_introns requires a sorted list, e.g. seq_100004 -> [93, 57] and seq_100004 -> [93, 62]
        # and retain only unique enries: seq_100009 -> 14.0 ... exon lengths: [426, 249, 369, 435, 0, 591, 390, 450, 246, 360, 357, 312, 321, 411] -> 0, because two times same position
        if len( self.position_list ) != len( set( self.position_list ) ):
            self.messages.append( '[ WARNING ] At least two intron positions are identical: {0}. This usually happens for short sequences. Although only one intron will be inserted, better take a closer look on the output!'.format( self.position_list ) )
        self.start_position_list = sorted( set( self.start_position_list ) )
        self.position_list = sorted( set( self.position_list ) )
        self.log_dict[ 'start_positions' ] = self.start_position_list
        self.log_dict[ 'insertion_positions' ] = self.position_list
        return



    def _determine_intermediate_positions( self, number_of_exons, start, remaining_len ):
        start_position_list = []
        intermediate_exon_len = int( round( remaining_len / number_of_exons ) )
        for i in range( 1, int(number_of_exons) ): #what if number of exons = 1 ???
            start_position = start + intermediate_exon_len * i
            start_position_list.append( start_position )
        return start_position_list, intermediate_exon_len



    def _get_insertion_positions(self, start_position, directions, max_search_range):
        position_list = []
        for i in range( max_search_range ):
            for sign in directions:
                position = start_position + i * sign
                two_n = self.dna_seq[ position : position + 2 ]
                if two_n == self.insertion_seq: # e.g., "GG"
                    position_list.append(position)
        position_list.reverse() # to use pop() method later
        return position_list
    
    
    
    def _is_cut_site_introduced( self, position, intron_seq, CSR_class ):
        if position - 10 < 0:
            downstream = self.dna_seq[ 0 : position + 1]
        else:
            downstream = self.dna_seq[ position - 10 : position + 1] #should end with 'G' for default settings
        upstream = intron_seq[ : 10 ]
        seq = downstream + upstream
        cut_site2index_list = CSR_class.get_cut_site_indices( dna_seq = seq )
        if cut_site2index_list:
            return cut_site2index_list
        # also consider end of intron seq
        downstream = intron_seq[ -10 : ]
        #if position + 11 > len(): # not required, since 'abcdef'[3:11] yields 'def'
        upstream = self.dna_seq[ position + 1 : position + 11 ] #should end with 'G' for default settings
        seq = downstream + upstream
        cut_site2index_list = CSR_class.get_cut_site_indices( dna_seq = seq )
        if cut_site2index_list:
            return cut_site2index_list        
        return False



    def insert_introns( self ):
        if not self.position_list:
            self.messages.append( '[ WARNING ] No positions to insert introns could be determined. Therefore, the output sequence does not contain any introns.' )
            l = [ self.dna_seq ]
        else:
            l = []
            if len( self.position_list ) == 1:
                position = self.position_list[ 0 ]
                l += [ self.dna_seq[ : position + 1 ], self.intron_seq, self.dna_seq[ position + 1 : ] ]
            else:
                for i, position in enumerate( self.position_list ):
                    if i == 0: # insert first intron -> take complete sequence from the beginning
                        l += [ self.dna_seq[ 0 : position + 1 ], self.intron_seq ]
                        continue
                    previous_position = self.position_list[ i - 1 ]
                    if i == len( self.position_list ) - 1: # last intron -> complete sequence with rest
                        l += [ self.dna_seq[ previous_position + 1 : position + 1 ], self.intron_seq, self.dna_seq[ position + 1 :  ] ]
                    else:
                        l += [ self.dna_seq[ previous_position + 1 : position + 1 ], self.intron_seq ]
        self.dna_seq_new = ''.join( l )
        return



    def check_intron_inserted_seq( self, CO_class, codon2aa_oneletter ):
        s = ''.join( self.dna_seq_new.split( self.intron_seq ) ).strip()
        return_value = True
        if CO_class.translate( self.dna_seq, codon2aa_oneletter ) != CO_class.translate( s, codon2aa_oneletter ):
            self.messages.append( '[ ERROR ] final DNA seq does not match input AA seq: {0}\n{1}\n{2}\n{3}\n{4}'.format(
                CO_class.translate( self.dna_seq, codon2aa_oneletter ), CO_class.translate( s, codon2aa_oneletter ), self.position_list, self.dna_seq, self.dna_seq_new
                    )
                )
            return_value = False
        if len( s ) % 3 != 0:
            self.messages.append( '[ ERROR ] final DNA seq length is not a multiplier of 3. {0}'.format(len(s)) )
            return_value = False
        if len( s ) == 0:
            self.messages.append( '[ ERROR ] final DNA seq length is 0. {0}'.format(len(s)) )
            return_value = False
        TMP = Tables( '' )
        if not set(s).issubset( TMP.allowed_characters_nt ):
            self.messages.append( '[ ERROR ]: final DNA seq contains invalid nt-characters. {0}'.format(sorted(set(s))) )
            return_value = False
        validation = set()
        for i, c in enumerate( s ):
            codon = False
            if i == 0:
                continue
            if i % 3 == 0:
                codon = s[ i - 3 : i ]
            if i == len( s ) - 1:
                codon = s[ i - 2 : i + 1 ]
            if codon:
                validation.add( TMP.check_codon( codon ) )
        if validation != set( [ True ] ):
            return_value = False
            self.messages.append( '[ ERROR ] final DNA seq contains invalid codons. {0}'.format( validation ) )
        return return_value



class MakeGenbank():
    
    def convert_position_list_into_start_end(self, position_list):
        assert len(position_list) > 1
        new_position_list = [ [ position_list[0] ] ]
        for i, pos in enumerate(position_list[1:]):
            if position_list[i] + 1 == pos:
                new_position_list[ -1 ].append(pos)
            else:
                new_position_list.append( [ pos ] )
        # check
        for l in new_position_list:
            assert sorted(l) == l
        # only retain start and end
        join_list = []
        for l in new_position_list:
            join_list.append( ( l[0] + 1, l[-1] + 1 ) ) # +1, because python-index-coordinates into gb-coordinates
        return join_list
    
    def convert_start_end_into_join(self, start_end_list):
        if not start_end_list:
            s = ''
        elif len(start_end_list) == 1:
            start, end = start_end_list[0]
            s = '{0}..{1}'.format(start, end)
        else:
            l = []
            for start, end in start_end_list:
                l.append('{0}..{1}'.format(start, end))
            s = 'join({0})'.format( ','.join(l) )
        return s
    
    def get_feature_string(self, position, label, feature_type = 'CDS', translation=''):
        if not position:
            s = ''
        elif not translation:
            s = '''     {feature_type}{position}
                     /label="{label}"'''
            s = s.format(feature_type = '{0:16}'.format(feature_type), position=position, label=label)
        else:
            s = '''     {feature_type}{position}
                     /label="{label}"
                     /translation="{translation}"'''
            s = s.format(feature_type = '{0:16}'.format(feature_type), position=position, label=label, translation=translation)
        return s

    def get_translation_string(self, aa_seq):
        breaks_list = [ 44 + 58 * i for i in range(100) ]
        translation = ''
        for i, aa in enumerate(aa_seq):
            if i in breaks_list:
                translation += os.linesep + ' '*21
            translation += aa
        return translation
    
    def _translate(self, dna_seq, codon2aa, CO_class):
        cds = SeqRecord( Seq(dna_seq, IUPAC.unambiguous_dna) )
        translation = CO_class.translate(dna_seq, codon2aa)
        assert cds.translate().seq == translation
        return translation

    def generate_gb_string(self, aa_seq, seqlist_annotated, codon2aa, CO_class, fasta_name = 'unnamed_input_seq', intron_name_seq_list = None):
        dna_seq_fine_tuned = ''.join([ n for n, annotation in seqlist_annotated ])
        if intron_name_seq_list == None:
            intron_name_seq_list = [ ('intron', '') ]
        start_list, linker_start_list, cds_list, intron_list, linker_end_list, stop_list = [], [], [], [], [], []
        cds_translation, linker_start_translation, linker_end_translation = '', '', ''
        for i, (n, annotation) in enumerate(seqlist_annotated):
            if annotation == 'Start':
                start_list.append(i)
            elif annotation == "5'-Linker":
                linker_start_list.append(i)
                linker_start_translation += n
            elif annotation == 'CDS':
                cds_list.append(i)
                cds_translation += n
            elif annotation == 'intron':
                intron_list.append(i)
            elif annotation == "3'-Linker":
                linker_end_list.append(i)
                linker_end_translation += n
            elif annotation == 'Stop':
                stop_list.append(i)
                
        if start_list:
            start_list = self.convert_position_list_into_start_end(start_list)
        if linker_start_list:
            linker_start_list = self.convert_position_list_into_start_end(linker_start_list)
            linker_start_translation = self.get_translation_string(
                self._translate(linker_start_translation, codon2aa=codon2aa, CO_class=CO_class))
        if cds_list:
            cds_list = self.convert_position_list_into_start_end(cds_list)
            cds_translation = self.get_translation_string(
                self._translate(cds_translation, codon2aa=codon2aa, CO_class=CO_class))
        if intron_list:
            intron_list = self.convert_position_list_into_start_end(intron_list)
        if linker_end_list: 
            linker_end_list = self.convert_position_list_into_start_end(linker_end_list)
            linker_end_translation = self.get_translation_string(
                self._translate(linker_end_translation, codon2aa=codon2aa, CO_class=CO_class))
        if stop_list:
            stop_list = self.convert_position_list_into_start_end(stop_list)
            
        # get intron names
        if intron_list:
            intron_name_list = [ intron_name_seq_list[0][0] for _ in intron_list ]
            if len(intron_name_seq_list) == 2:
                intron_name_list[-1] = intron_name_seq_list[1][0]

        # prepare dna string
        origin = ''
        for i, nt in enumerate( dna_seq_fine_tuned ):
            if i % 60 == 0:
                if i == 0:
                    tmp = ''
                else:
                    tmp = os.linesep
                origin += tmp + '{0:9}'.format(i + 1)
            if i % 10 == 0:
                origin += ' '
            origin += nt
        
        # gb string
        s = '''LOCUS       {scaffold}           {length} bp    DNA     PLN       {date}
DEFINITION  ChlamyIntronserter optimized cDNA with introns
ACCESSION   None
VERSION     None.1
KEYWORDS    .
SOURCE      Generated by the ChlamyIntronserter tool
  ORGANISM  unknown
REFERENCE   1  (bases 1 to {length})
  AUTHORS   Baier,T., Wichmann,J., Kruse,O. and Lauersen,K.J.,
  TITLE     Intron-containing algal transgenes mediate efficient
            recombinant gene expression in the green microalga
            Chlamydomonas reinhardtii
  JOURNAL   Nucleic Acids Research. doi:10.1093/nar/gky532. (2018)
REFERENCE   2  (bases 1 to {length})
  AUTHORS   Jaeger,D. and Lauersen,K.J.
  TITLE     tba
  JOURNAL   tba (2018)
FEATURES             Location/Qualifiers
     CDS             {cds_list}
                     /label="{fasta_name}__CDS"
                     /codon_start=1
                     /product="{fasta_name}"
                     /translation="{translation}"
{intron_string}
{linker_start}
{linker_end}
{start}
{stop}
ORIGIN
{origin}
//
'''
        genbank_string = s.format(
            translation = cds_translation,
            origin = origin,
            fasta_name = fasta_name,
            length = len(dna_seq_fine_tuned),
            cds_list = self.convert_start_end_into_join(cds_list),
            scaffold = 'CIoptDNA',
            date = '{0}'.format( datetime.date.today() ),
            intron_string = os.linesep.join(
                [ self.get_feature_string(position=self.convert_start_end_into_join([intron]), label=intron_name_list[i], feature_type='intron')
                  for i, intron in enumerate(intron_list) ]
            ),
            linker_start = self.get_feature_string(
                position=self.convert_start_end_into_join(linker_start_list), label='5-linker', translation=linker_start_translation),
            linker_end = self.get_feature_string(
                position=self.convert_start_end_into_join(linker_end_list), label='3-linker', translation=linker_end_translation),
            start = self.get_feature_string(
                position=self.convert_start_end_into_join(start_list), label='start', translation='M'),
            stop = self.get_feature_string(
                position=self.convert_start_end_into_join(stop_list), label='stop', translation='*'),
            )
        genbank_string = os.linesep.join([s for s in genbank_string.splitlines() if s]) #remove empty lines from string
        return genbank_string
    
    
    
    def check_gb(self, aa_seq, gb_string):
        with tempfile.TemporaryFile(mode='w+') as fp: # + for mode open a disk file for updating (reading and writing) -> https://docs.python.org/3/library/functions.html#open
            fp.write(gb_string)
            fp.seek(0)
            seq_record = SeqIO.read(fp, "genbank")
            translation_list = []
            start_end_list = []
            for feature in seq_record.features:
                if feature.type == 'CDS':
                    if 'translation' in feature.qualifiers:
                        assert len(feature.qualifiers['translation']) == 1
                        translation_list.append( (feature.location.start, feature.qualifiers['translation'][0]) )
                    for exon in feature.location.parts:
                        start_end_list.append( (exon.start, exon.end) )
            start_end_list.sort()
            cds = SeqRecord( Seq('', IUPAC.unambiguous_dna) )
            for start, end in start_end_list:
                cds += seq_record.seq[ start : end ]
        translation_list.sort()
        translation = ''.join([ translation for start, translation in translation_list ])
        if cds.translate().seq == aa_seq and translation == aa_seq:
            return True
        else:
            return False



class PrepareLog():

    def cut_site_removal_log( self, log_dict ):
        output_string = ''
        if not log_dict:
            output_string = '[ INFO ] No cut site(s) in sequence detected.'
        for cut_site in sorted( log_dict ):
            output_string += '<h3>Detected cut site: {0}</h3>'.format( cut_site )
            for position in log_dict[ cut_site ]:
                start_bold_font = position - log_dict[ cut_site ][ position ][ 'index_2_start_from' ]
                output_string += '<p>... detected at position {3} in cDNA sequence, region with cut site: {0}<b>{1}</b>{2}; all synonymous codon replacement possibilities:</p>'.format(
                    log_dict[ cut_site ][ position ][ 'region_with_cut_site' ][ : start_bold_font ],
                    log_dict[ cut_site ][ position ][ 'region_with_cut_site' ][start_bold_font : start_bold_font+len(cut_site) ],
                    log_dict[ cut_site ][ position ][ 'region_with_cut_site' ][ start_bold_font+len(cut_site) : ],
                    position )

                output_string += '<table> <tr class="tr_alternate_color">'
                output_string += '<th>{0}</th> <th>{1}</th> <th>{2}</th> <th>{3}</th> </tr>'.format( 'original codon&nbsp;&nbsp;&nbsp;', 'alternative codon&nbsp;&nbsp;&nbsp;', 'frequency&nbsp;&nbsp;&nbsp;', 'cut site removed?' )
                for freq, possible_codon, old_codon, cut_site_removed in log_dict[ cut_site ][ position ][ 'freq_altcodon_orgcodon_remove_list' ]:
                    if possible_codon == log_dict[ cut_site ][ position ][ 'new_codon' ]:
                        output_string += '<tr class="tr_alternate_color"> <td><b>{0}</b></td> <td><b>{1}</b></td> <td><b>{2}</b></td> <td><b>{3}</b></td> </tr>'.format( old_codon, possible_codon, freq, cut_site_removed )
                    else:
                        if cut_site_removed:
                            output_string += '<tr class="tr_alternate_color"> <td>{0}</td> <td>{1}</td> <td>{2}</td> <td>{3}</td> </tr>'.format( old_codon, possible_codon, freq, cut_site_removed )
                        else:
                            output_string += '<tr class="tr_alternate_color"> <td style="color:Grey">{0}</td> <td style="color:Grey">{1}</td> <td style="color:Grey">{2}</td> <td style="color:Grey">{3}</td> </tr>'.format( old_codon, possible_codon, freq, cut_site_removed )
                output_string += '</table>'
                output_string += '<p>region without cut site: {0}<b>{1}</b>{2}</p>'.format(
                    log_dict[ cut_site ][ position ][ 'region_without_cut_site' ][ : start_bold_font ],
                    log_dict[ cut_site ][ position ][ 'region_without_cut_site' ][ start_bold_font : start_bold_font+len(cut_site) ],
                    log_dict[ cut_site ][ position ][ 'region_without_cut_site' ][ start_bold_font+len(cut_site) : ] )
        return output_string

    def intron_insertion_log( self, log_dict ):
        output_string = '<p><ul>'
        if 'seq_len' in log_dict: # automatic determination
            output_string += '<li>Automatic determination of start positions was performed, because start positions were not explicitly given.</li>'
            if log_dict[ 'initial_exon_determination' ] == None:
                if len( log_dict[ 'start_positions' ] ) == 1:
                    output_string += '<li><b>Exception for short sequences:</b> Only the end-intron was inserted.</li>'
                else:
                    output_string += '<li><b>Exception for short sequences:</b> Only the start- and end-intron were inserted.</li>'
            else:
                min_exon_number, dist_1, max_exon_number, dist_2, number_of_exons = log_dict[ 'initial_exon_determination' ]
                output_string += '<li>The cDNA sequence has {0} nt and after insertion of a start-intron and an end-intron, the remaining length is {1} nt.</li>'.format(log_dict[ 'seq_len' ], log_dict[ 'seq_len_remaining' ])
                output_string += '<li>Dividing the sequence into {0} intermediate-<u>exons</u> (thus distributing {1} intermediate-<u>introns</u> into the remaining sequence) results in a difference to the parameter <span style="font-family:\'Courier New\'">intermediate</span> by {2} nt.</li>'.format(min_exon_number, min_exon_number-1, dist_1)
                output_string += '<li>Dividing the sequence into {0} intermediate-<u>exons</u> (thus distributing {1} intermediate-<u>introns</u> into the remaining sequence) results in a difference to the parameter <span style="font-family:\'Courier New\'">intermediate</span> by {2} nt.</li>'.format(max_exon_number, max_exon_number-1, dist_2)
                output_string += '<li>In order to minimize the distance to the parameter <span style="font-family:\'Courier New\'">intermediate</span>, the initial number of intermediate-<u>exons</u> is defined as {0}.</li>'.format(number_of_exons)
                target_possible, intermediate_exon_len, final_number_of_exons = log_dict[ 'final_exon_determination' ]
                if target_possible:
                    output_string += '<li>With the initially defined number of intermediate-<u>exons</u>, the lengths of intermediate-exons are consistently below the parameter <span style="font-family:\'Courier New\'">max exon length</span>. Therefore, {0} intermediate-<u>introns</u> were inserted into the sequence.</li>'.format(number_of_exons - 1)
                else:
                    output_string += '<li>With the initially defined number of {0} intermediate-<u>exons</u>, the maximal length of intermediate-exons exceeded the parameter <span style="font-family:\'Courier New\'">max exon length</span>. Choosing {1} intermediate-exons instead, the lengths of intermediate-exons were found to remain consistently below <span style="font-family:\'Courier New\'">max exon length</span>. Therefore, {2} intermediate-<u>introns</u> were inserted into the sequence.</li>'.format(number_of_exons, final_number_of_exons, final_number_of_exons - 1)
                target_intermediate_exon_len = log_dict[ 'target_intermediate_exon_len' ]
                output_string += '<li>With the finally defined number of {0} intermediate-<u>introns</u>, a target intermediate-<u>exon</u> length of {1} nt had to be used (instead of {2} nt as requested by the parameter <span style="font-family:\'Courier New\'">intermediate</span>).</li>'.format(final_number_of_exons-1, intermediate_exon_len, target_intermediate_exon_len)
        else: # manual input
            output_string += '<li>The start positions for the insertion of introns were manually given as input. Therefore, an automatic determination of start positions was not required and thus not performed.</li>'
        output_string += '<li>Final "start positions" for the insertion of introns were: {0}.</li>'.format( log_dict[ 'start_positions' ] )
        output_string += '<li>Actual insertion positions of the introns were: {0}.</li>'.format( log_dict[ 'insertion_positions' ] )
        output_string += '</ul></p>'
        return output_string


class Plotting():

    def plot_norm_codon_freq(self, dna_seq, aa2codon_freq_list):
        # convert dna_seq to codon_list:
        codon = ''
        codon_list = []
        for i, nt in enumerate( dna_seq ):
            codon += nt
            if (i + 1) % 3 == 0:
                codon_list.append( codon )
                codon = ''
        # normalize frequencies:
        codon2norm_freq = {}
        codon2aa = {}
        for aa, l in aa2codon_freq_list.items():
            if len( l ) == 1:
                codon, freq = l[ 0 ]
                codon2norm_freq[ codon ] = 100
                continue
            freq_list = [ freq for codon, freq in l ]
            max_freq = max( freq_list )
            min_freq = min( freq_list )
            for i, (codon, freq) in enumerate( l ):
                if max_freq - min_freq == 0: # both have to be 0.5
                    assert len(l) == 2
                    codon2norm_freq[ codon ] = 100 # replacing a 0.5 codon by a 0.5 codon yields a norm freq of 100
                else:
                    codon2norm_freq[ codon ] = (freq - min_freq) / (max_freq - min_freq) * 100
                codon2aa[ codon ] = aa
        # prepare lists
        x, y = [], []
        for i, codon in enumerate( codon_list ):
            x.append( i )
            y.append( codon2norm_freq[ codon ] )
        plt.xlabel('codon in cDNA sequence [position in aa coordinates]')
        plt.ylabel('normalized codon frequency [%]')
        plt.title( 'Normalized Codon Frequencies of the cDNA sequence\n(values < 100% are due to restriction digest removal by synonymous codon replacement)\n\n' )
        plt.plot( x, y, '-b' ) # add blue lines
        for i, codon in enumerate( codon_list ):
            if codon2norm_freq[ codon ] != 100:
                plt.text( i, codon2norm_freq[ codon ] - 5, '{0}'.format( codon ), horizontalalignment='center', size='x-small' )
                for codon_2, freq in aa2codon_freq_list[ codon2aa[ codon ] ]:
                    plt.plot( i, codon2norm_freq[ codon_2 ], 'or' )

        plt.legend( [ 'norm. codon freq.', 'available codons' ], bbox_to_anchor = (0.0, 1.02, 1., 0.), loc=4, ncol=2, borderaxespad=0. )
        fig = plt.gcf()
        fig.set_size_inches(12, 4)

        figfile = BytesIO()
        plt.savefig(figfile, format='png', bbox_inches='tight') # tight removes any white space
        #plt.savefig("norm_freq.svg")
        figfile.seek(0)  # rewind to beginning of file
        figdata_png = base64.b64encode(figfile.getvalue())
        figdata_png = figdata_png.decode('utf8')

        plt.close()
        self.normfreq_png = figdata_png
        return figdata_png

    def plot_gene_architecture(self, dna_seq_list, cut_sites):
        cut_site_start, cut_site_end = cut_sites
        y_pos = [ 1, 1 ]
        lwd_ex, lwd_in = 6.0, 3.0
        #col_ex, col_in = '#B7DEE8', '#FAC090'
        col_ex, col_in = '#9ABA59', '#708BCD'

        length_list = [ len(ex_in) for ex_in in dna_seq_list ]
        # first draw introns
        for i, ex_in in enumerate( dna_seq_list ):
            if i % 2 == 0: # exons
                continue
            else: # introns
                intron = ex_in
                intron_start_position = sum(length_list[:i])
                intron_end_position = sum(length_list[:i]) + len( intron ) - 1
                plt.plot( [ intron_start_position, intron_end_position ], y_pos, linestyle='-', color=col_in, linewidth=lwd_in )
        # second draw exons
        for i, ex_in in enumerate( dna_seq_list ):
            if i % 2 == 0: # exons
                exon = ex_in
                exon_len = len( exon )
                exon_start_position = sum(length_list[:i]) + 1
                exon_end_position = sum(length_list[:i]) + exon_len
                if i == 0: # first exon
                    exon_start_position = 1 + len( cut_site_start )
                    exon_len -= len( cut_site_start )
                elif i == len( dna_seq_list ) - 1: # last exon
                    exon_end_position = sum(length_list[:i]) + len( exon ) - len(cut_site_end)
                    exon_len -= len(cut_site_end)
                plt.plot( [ exon_start_position, exon_end_position ], y_pos, color=col_ex, linestyle='-', linewidth=lwd_ex )
                plt.text( exon_start_position + exon_len / 2., 1.005, '{0}'.format( exon_len ), horizontalalignment='center', size='small' )
            else: # introns
                continue

        plt.xlabel('position in optimized DNA sequence [nt]')
        plt.ylabel('')
        plt.title( 'Exon-Intron structure of the optimized DNA sequence' )
        blue_line = mlines.Line2D( [], [], color = col_ex, label = 'exon')
        green_line = mlines.Line2D( [], [], color = col_in, label = 'intron')
        plt.legend( handles = [ blue_line, green_line ] )
        ax = plt.gca()
        #ax.yaxis.set_visible(False)
        #ax.spines['left'].set_color('white')
        ax.yaxis.label.set_color('white')
        ax.tick_params(axis='y', colors='white')
        fig = plt.gcf()
        fig.set_size_inches(12, 2)
        figfile = BytesIO()
        plt.savefig(figfile, format='png', bbox_inches='tight') # tight removes any white space
        #plt.savefig("in_ex.svg")
        figfile.seek(0)  # rewind to beginning of file
        figdata_png = base64.b64encode(figfile.getvalue())
        figdata_png = figdata_png.decode('utf8')
        plt.close()
        self.exInStructure_png = figdata_png
        return figdata_png
