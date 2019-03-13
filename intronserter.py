#!/usr/bin/env python3

#BSD 3-Clause License
#
#Copyright (c) 2019, Daniel Jaeger
#All rights reserved.
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are met:
#
#* Redistributions of source code must retain the above copyright notice, this
#list of conditions and the following disclaimer.
#
#* Redistributions in binary form must reproduce the above copyright notice,
#this list of conditions and the following disclaimer in the documentation
#and/or other materials provided with the distribution.
#
#* Neither the name of the copyright holder nor the names of its
#contributors may be used to endorse or promote products derived from
#this software without specific prior written permission.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
#CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
#OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import os
import sys
import collections
import string
import argparse
import base64
import traceback
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import class_library

class MessageContainer():
    def __init__(self):
        self.messages = collections.OrderedDict()
        self.messages[ 'global' ] = []
        return

def parse_input( ArgsClass, MessageContainer ):
    TMP = class_library.Tables( '' )
    intron_name_seq_list = []

    # update 13.03.2019: SeqIO.parse does not iterate when FASTA file contains no ">name" as first line, but only seq ...
    default_name = 'unnamed_input_seq'
    invalid_fasta = False
    with open(ArgsClass.aa_fasta_file, 'rU') as fin:
        for i, line in enumerate(fin):
            if i == 0:
                if not line.startswith(">"):
                    invalid_fasta = True
    if invalid_fasta:
        with open(ArgsClass.aa_fasta_file, 'rU') as fin:
            s = fin.read()
        with open(ArgsClass.aa_fasta_file, 'w') as fout:
            print('>{0}'.format(default_name), file=fout)
            print(s, file=fout)

    # validate the parsed fasta seq by comparison against Biopython
    with open(ArgsClass.aa_fasta_file, 'rU') as fin:
        aa_seq_dict = TMP.parse_fasta( fin.read(), default_name=default_name )
        fin.seek(0)
        for (name, seq), record in zip(aa_seq_dict.items(), SeqIO.parse(fin, "fasta")):
            #record.name == name is FALSE when the header contains any white spaces - white spaces are removed by Biopython
            if not record.seq.upper() == seq:
                MessageContainer.messages['global'].append('[ ERROR ] Parsing the FASTA AA seq of {0} was unsuccessful. VALIDATE output!'.format(name))

    # validate FASTA AA input for only valid characters
    allowed_characters = set(IUPAC.protein.letters + '*')
    assert sorted([ TMP.AA_lookup_dict[k]['1-letter'] for k in TMP.AA_lookup_dict ]) == sorted(allowed_characters)
    allowed_characters_DNA = IUPAC.unambiguous_dna.letters
    for name, seq in aa_seq_dict.items():
        MessageContainer.messages[name] = []
        if ArgsClass.only_insert_introns:
            if not set(seq).issubset(allowed_characters_DNA):
                MessageContainer.messages[name].append('[ ERROR ] You specified --only_insert_introns, but your FASTA sequence {0} contains invalid characters. Allowed characters for this option are {1}'.format(name, sorted(allowed_characters_DNA)))
                aa_seq_dict[name] = ''
        else:
            if not set(seq).issubset(allowed_characters):
                MessageContainer.messages[name].append('[ ERROR ] Your FASTA sequence {0} contains invalid characters. Allowed characters are {1}'.format(name, sorted(allowed_characters)))
                aa_seq_dict[name] = ''

    # check if aa- or cDNA-seq:
    allowed_characters = set('ATCG')
    for name, seq in aa_seq_dict.items():
        if seq and set(seq).issubset(allowed_characters):
            if len(seq) % 3 != 0:
                MessageContainer.messages[name].append( '[ ERROR ] Your input sequence "{0}" seems to be a cDNA sequence, but its length is not a multiplier of 3! Please re-submit a valid sequence, preferably as amino acids.'.format(name) )
                aa_seq_dict[ name ] = ''
            else:
                if ArgsClass.only_insert_introns:
                    MessageContainer.messages[name].append( '[ INFO ] Your input sequence "{0}" seems to be a cDNA sequence and the option --only_insert_introns was specified. Only introns will be inserted, nothing else.'.format(name) )
                else:
                    MessageContainer.messages[name].append( '[ INFO ] Your input sequence "{0}" seems to be a cDNA sequence, because it contains only the characters A, T, C and G. The sequence will be translated to its amino acid counterpart first.'.format(name) )
                    aa_seq_dict[ name ] = class_library.CodonOptimizer( '' ).translate( seq, TMP.codon2aa_oneletter )

    if ArgsClass.custom_codon_usage_table_file:
        with open(ArgsClass.custom_codon_usage_table_file) as fin:
            codon_usage_table = fin.read()
            try:
                TMP.import_codon_table( codon_table_string = codon_usage_table )
                TMP.convert_codon_counts_2_freq()
            except:
                tb = traceback.format_exc()
                MessageContainer.messages['global'].append('[ ERROR ] Invalid codon table. Using default=Kazusa for C. reinhardtii instead... error message was: {0}'.format(tb))
                codon_usage_table = TMP._kazusa_codon_table()
    else:
        if ArgsClass.codon_usage_table_id == 'kazusa':
            codon_usage_table = TMP._kazusa_codon_table() # get internal table
        else:
            codon_usage_table = TMP._hivecut_codon_table() # get internal table

    cut_site_list = [ cut_site for cut_site in ArgsClass.cut_sites.split(',') if cut_site ]
    if 'custom' in cut_site_list:
        cut_site_list = [ cut_site for cut_site in cut_site_list if cut_site != 'custom' ] # remove entry 'custom'
        if ArgsClass.custom_cut_sites:
            cut_site_list += [ cut_site for cut_site in ArgsClass.custom_cut_sites.split(',') if cut_site ] # process cut site fasta input
    allowed_len = set([6, 8])
    for cut_site in cut_site_list:
        if len(cut_site) not in allowed_len:
            MessageContainer.messages['global'].append('[ ERROR ] The list of cut sites to avoid contains at least one sequence with a length unequal to either 6 or 8. This/these sequences are ignored.')
            break
    cut_site_list = [cut_site for cut_site in cut_site_list if len(cut_site) == 6 or len(cut_site) == 8]

    intron_name_seq_list.append( ('intron', ArgsClass.intron_seq.lower()) )

    if ArgsClass.intron_lastdifferent:
        if ArgsClass.intron_lastdifferent_seq:
            intron_name_seq_list.append( ( 'last_intron', ArgsClass.intron_lastdifferent_seq.lower() ) )
        else:
            MessageContainer.messages['global'].append('[ ERROR ] You specified the parameter "--intron_lastdifferent", but did not specify the parameter "--intron_lastdifferent_seq". Option is ignored and the last intron will NOT be substituted. ')

    if ArgsClass.supersede_intron_insert:
        if ArgsClass.manual_intron_positions:
            start_position_list = [ int(position.strip()) for position in ArgsClass.manual_intron_positions.split(',') if position ]
        else:
            MessageContainer.messages['global'].append('[ ERROR ] You specified the parameter "--supersede_intron_insert", but did not specify the parameter "--manual_intron_positions". Option is ignored and AUTOMATIC intron insertion performed. ')
            start_position_list = None
    else:
        start_position_list = None

    # cut site
    cut_site_start = ArgsClass.cut_site_start
    if cut_site_start == 'None':
        cut_site_start = ''
    elif cut_site_start == 'custom':
        cut_site_start = ArgsClass.custom_cut_site_start

    cut_site_end = ArgsClass.cut_site_end
    if cut_site_end == 'None':
        cut_site_end = ''
    elif cut_site_end == 'custom':
        cut_site_end = ArgsClass.custom_cut_site_end

    kwargs = {}
    kwargs[ 'aa_seq_dict' ] = aa_seq_dict
    kwargs[ 'codon_table' ] = codon_usage_table
    kwargs[ 'cut_site_list' ] = cut_site_list
    kwargs[ 'intron_seq' ] = ArgsClass.intron_seq.lower()
    kwargs[ 'intron_name_seq_list' ] = intron_name_seq_list
    kwargs[ 'insertion_seq' ] = ArgsClass.nucleotide_pair
    kwargs[ 'start' ] = ArgsClass.start
    kwargs[ 'intermediate' ] = ArgsClass.target
    kwargs[ 'end' ] = ArgsClass.end
    kwargs[ 'max_exon_length' ] = ArgsClass.max
    kwargs[ 'start_position_list' ] = start_position_list
    
    kwargs[ 'cut_site_start' ] = cut_site_start.upper()
    kwargs[ 'cut_site_end' ] = cut_site_end.upper()
    kwargs[ 'linker_start' ] = ArgsClass.linker_start.upper()
    kwargs[ 'linker_end' ] = ArgsClass.linker_end.upper()
    kwargs[ 'insert_start_codon' ]  = True if ArgsClass.insert_start_codon else False
    kwargs[ 'insert_stop_codon' ]   = True if ArgsClass.insert_stop_codon else False
    kwargs[ 'remove_start_codon' ]  = True if ArgsClass.remove_start_codon else False
    kwargs[ 'remove_stop_codon' ]   = True if ArgsClass.remove_stop_codon else False

    kwargs[ 'only_insert_introns' ] = True if ArgsClass.only_insert_introns else False

    kwargs[ 'output_dict' ] = None
    return kwargs



def process_input( kwargs, MessageContainer ):
    output_dict = collections.OrderedDict()
    remove_punctuation_map = dict((ord(char), None) for char in string.punctuation)

    # check if the intron sequences are free from cut sites
    if not kwargs[ 'only_insert_introns' ]:
        CSR_class = class_library.CutSiteRemover( dna_seq = '', cut_site_list = kwargs[ 'cut_site_list' ], codon2aa = {}, aa2codon_freq_list = {} )
        for name, seq in kwargs[ 'intron_name_seq_list' ]:
            seq = kwargs[ 'insertion_seq' ][0] + seq + kwargs[ 'insertion_seq' ][1]
            cut_site2index_list = CSR_class.get_cut_site_indices( dna_seq = seq )
            if cut_site2index_list:
                l = []
                for cut_site, index_found_at_list in cut_site2index_list.items():
                    l.append( 'cut site "{0}" at the position(s): {1}'.format( cut_site, index_found_at_list ) )
                MessageContainer.messages['global'].append( '[ ERROR ] The following cut sites are part of the "{1}" intron sequence: {0}. (a) Intron insertion might NOT be possible and (b) these cut sites are NOT removed from the optimized sequence.'.format( l, name ) )

    # fine tune sequence
    if not kwargs[ 'only_insert_introns' ]:
        start_codon, stop_codon = '', ''
        if kwargs[ 'insert_start_codon' ]:
            start_codon = 'M'
        if kwargs[ 'insert_stop_codon' ]:
            stop_codon = '*'
            if kwargs[ 'linker_end' ]:
                MessageContainer.messages['global'].append( '[ INFO ] Although you requested to insert a * stop codon, the * stop codon was NOT inserted, because you also requested a 3\'-linker. Inserting the * stop codon would have resulted in translation termination, counteracting your intended protein fusion as indicated by the 3\'-linker insertion request.' )
                kwargs[ 'insert_stop_codon' ] = False
                stop_codon = ''
        for name, aa_seq in kwargs[ 'aa_seq_dict' ].items():
            if not aa_seq:
                continue
            aa_seq = aa_seq.upper()
            # remove start codon if requested
            for  z in range(1000):
                if kwargs[ 'remove_start_codon' ] and aa_seq.startswith('M'):
                    aa_seq = aa_seq[1:]
                else:
                    break
            # remove stop codon if requested or if 3'-linker = linker_end is given
            for  z in range(1000):
                if (kwargs[ 'remove_stop_codon' ] or kwargs[ 'linker_end' ]) and aa_seq.endswith('*'):
                    aa_seq = aa_seq[:-1]
                    if kwargs[ 'linker_end' ] and not kwargs[ 'remove_stop_codon' ]:
                        MessageContainer.messages['global'].append( '[ INFO ] Although you did not request to remove the native * stop codon, the * stop codon was automatically removed, because you requested to insert a 3\'-linker. Not removing the * stop codon would have resulted in translation termination, counteracting your intended protein fusion as indicated by the 3\'-linker insertion request.' )
                else:
                    break
            # insert start codon, linker_start, linker end, stop_codon
            kwargs[ 'aa_seq_dict' ][ name ] = start_codon + kwargs[ 'linker_start' ] + aa_seq + kwargs[ 'linker_end' ] + stop_codon
        
    for j, (name, aa_seq) in enumerate( kwargs[ 'aa_seq_dict' ].items() ):
        if not aa_seq:
            continue
        # initialize:
        if kwargs[ 'only_insert_introns' ]:
            DB_class = class_library.Tables( str(Seq(aa_seq, IUPAC.unambiguous_dna).translate()).upper() )
        else:
            DB_class = class_library.Tables( aa_seq.upper() )
        output_dict[ name ] = { 'cDNA_seq_plus_i' : None, 'genbank_string' : None, 'name' : '>{0}'.format(name) }

        # codon optimize:
        if not kwargs[ 'only_insert_introns' ]:
            DB_class.import_codon_table( codon_table_string = kwargs[ 'codon_table' ] )
            DB_class.convert_codon_counts_2_freq()
            CO_class = class_library.CodonOptimizer( aa_seq = DB_class.aa_seq )
            CO_class.reverse_translate( DB_class.aa_oneletter_2_mostfreq_codon_dict )
            DB_class.cDNA_seq = CO_class.dna_seq
        else:
            DB_class.import_codon_table( codon_table_string = kwargs[ 'codon_table' ] )
            DB_class.convert_codon_counts_2_freq()
            CO_class = class_library.CodonOptimizer( aa_seq = '' )

        # cut site removal:
        if not kwargs[ 'only_insert_introns' ]:
            CSR_class = class_library.CutSiteRemover( dna_seq = DB_class.cDNA_seq, cut_site_list = kwargs[ 'cut_site_list' ], codon2aa = DB_class.codon2aa, aa2codon_freq_list = DB_class.aa2codon_freq_list )
            DB_class.cDNA_seq_cleaned = CSR_class.main( iter_max = 1000 )
            MessageContainer.messages[name] += CSR_class.messages
        else:
            CSR_class = class_library.CutSiteRemover( dna_seq = '', cut_site_list = [], codon2aa = DB_class.codon2aa, aa2codon_freq_list = DB_class.aa2codon_freq_list )
            DB_class.cDNA_seq_cleaned = aa_seq.upper() # aa_seq is in this case a cDNA seq
        
        # annotate the sequence using a mapping nucleotide->annotation
        seqlist = list( DB_class.cDNA_seq_cleaned )
        seqlist_annotated_beforeCDS, seqlist_annotated_afterCDS = [], []
        if kwargs['insert_start_codon']:
            for i in range(3):
                seqlist_annotated_beforeCDS.append( (seqlist.pop(0), 'Start') )
        if kwargs['linker_start']:
            for i in range(len(kwargs['linker_start'])*3):
                seqlist_annotated_beforeCDS.append( (seqlist.pop(0), "5'-Linker") )
        if kwargs['insert_stop_codon']:
            for i in range(3):
                seqlist_annotated_afterCDS.append( (seqlist.pop(), 'Stop') )
        if kwargs['linker_end']:
            for i in range(len(kwargs['linker_end'])*3):
                seqlist_annotated_afterCDS.append( (seqlist.pop(), "3'-Linker") )
        seqlist_annotated = seqlist_annotated_beforeCDS + [ (n, 'CDS') for n in seqlist ] + list(reversed(seqlist_annotated_afterCDS))
        assert len(seqlist_annotated) == len(DB_class.cDNA_seq_cleaned)
        assert ''.join([ n for n, a in seqlist_annotated ]) == DB_class.cDNA_seq_cleaned

        # intron insertion:
        if len(kwargs[ 'intron_name_seq_list' ]) == 2:
            intron_seq_2 = kwargs[ 'intron_name_seq_list' ][1][1]
        else:
            intron_seq_2 = ''
        II_class = class_library.IntronInserter(
            dna_seq = DB_class.cDNA_seq_cleaned,
            insertion_seq = kwargs[ 'insertion_seq' ],
            intron_seq = kwargs[ 'intron_seq' ],
            start = kwargs[ 'start' ],
            intermediate = kwargs[ 'intermediate' ],
            end = kwargs[ 'end' ],
            max = kwargs[ 'max_exon_length' ],
            start_position_list = kwargs[ 'start_position_list' ],
            end_intron_different_seq = intron_seq_2,
            CSR_class = CSR_class
            )
        II_class.determine_positions()
        II_class.insert_introns()
        MessageContainer.messages[name] += II_class.messages
        DB_class.cDNA_seq_plus_i = II_class.dna_seq_new
        
        # substitute the last intron to rbcS2 i2 if requested:
        dna_seq_list = []
        exon_list = DB_class.cDNA_seq_plus_i.split( II_class.intron_seq )
        last_intron_different = True if len( kwargs['intron_name_seq_list'] ) == 2 else False
        for i, exon in enumerate( exon_list ):
            dna_seq_list.append( exon )
            if i < len( exon_list ) - 1:
                if i == len( exon_list ) - 2 and last_intron_different:
                    intron_2_name, intron_2_seq = kwargs[ 'intron_name_seq_list' ][1]
                    dna_seq_list.append( intron_2_seq )
                else:
                    dna_seq_list.append( II_class.intron_seq )
        DB_class.cDNA_seq_plus_i = ''.join( dna_seq_list )
        
        # include introns into annotated seq
        for i, n in enumerate(DB_class.cDNA_seq_plus_i):
            if not n == seqlist_annotated[i][0]:
                seqlist_annotated.insert(i, ( n, 'intron' ))
        assert ''.join([ n for n, a in seqlist_annotated ]) == DB_class.cDNA_seq_plus_i

        # final check: no cut sites appeared due to insertion of introns
        if not kwargs[ 'only_insert_introns' ]:
            cut_site2index_list = CSR_class.get_cut_site_indices( dna_seq = DB_class.cDNA_seq_plus_i )
            if cut_site2index_list:
                l = []
                cut_site2enzyme = { v:k for k, v in DB_class.re_lookup_dict.items() }
                cut_site2enzyme.update( { CO_class.reverse_complement( v ) : k for k, v in DB_class.re_lookup_dict.items() } )
                for cut_site, index_found_at_list in cut_site2index_list.items():
                    l.append( 'cut site "{0}" ({1}) at the position(s): {2}'.format( cut_site, cut_site2enzyme[cut_site], index_found_at_list ) )
                MessageContainer.messages[name].append( '[ WARNING ] The following cut sites appeared due to the insertion of introns: {0}'.format( l ) )
            
        # final check 2: spliced, translated seq is identical to input seq
        exon_list = DB_class.cDNA_seq_plus_i.split( II_class.intron_seq )
        if last_intron_different:
            intron_2_name, intron_2_seq = kwargs[ 'intron_name_seq_list' ][1]
            last_two_exons = exon_list[-1].split(intron_2_seq)
            exon_list = exon_list[:-1] + last_two_exons
        coding_dna = ''.join(exon_list)
        if not kwargs[ 'only_insert_introns' ]:
            if not CO_class.translate(coding_dna, DB_class.codon2aa_oneletter) == aa_seq or not str(Seq(coding_dna, IUPAC.unambiguous_dna).translate()) == aa_seq:
                MessageContainer.messages[name].append( '[ ERROR ] The translation of the spliced optimized DNA sequence does not match the input AA sequence. MANUALLY VALIDATE OUTPUT!' )
        else:
            translated_intron_enriched_seq = CO_class.translate(coding_dna, DB_class.codon2aa_oneletter)
            translated_input_seq = CO_class.translate(aa_seq, DB_class.codon2aa_oneletter)
            translated_intron_enriched_seq_2 = str(Seq(coding_dna, IUPAC.unambiguous_dna).translate())
            translated_input_seq_2 = str(Seq(aa_seq, IUPAC.unambiguous_dna).translate())
            if not translated_intron_enriched_seq == translated_input_seq or not translated_intron_enriched_seq_2 == translated_input_seq_2:
                MessageContainer.messages[name].append( '[ ERROR ] The translation of the spliced optimized DNA sequence does not match the translation of the input DNA sequence. MANUALLY VALIDATE OUTPUT!' )


        # add cut sites:
        if kwargs[ 'only_insert_introns' ]:
            dna_seq_fine_tuned = DB_class.cDNA_seq_plus_i
            dna_seq_list_fine_tuned = list(dna_seq_list)
        else:
            dna_seq_fine_tuned = kwargs[ 'cut_site_start' ] + DB_class.cDNA_seq_plus_i  + kwargs[ 'cut_site_end' ]
            dna_seq_list_fine_tuned = list(dna_seq_list)
            dna_seq_list_fine_tuned[ 0 ] = kwargs[ 'cut_site_start' ] + dna_seq_list_fine_tuned[ 0 ]
            dna_seq_list_fine_tuned[ -1 ] = kwargs[ 'cut_site_end' ] + dna_seq_list_fine_tuned[ -1 ]
        
            # include cut sites into annotated seq
            if kwargs[ 'cut_site_start' ]:
                seqlist_annotated = [ (n, 'cut_site') for n in kwargs[ 'cut_site_start' ] ] + seqlist_annotated
            if kwargs[ 'cut_site_end' ]:
                seqlist_annotated = seqlist_annotated + [ (n, 'cut_site') for n in kwargs[ 'cut_site_end' ] ]
            assert ''.join([ n for n, a in seqlist_annotated ]) == dna_seq_fine_tuned
            
        # create genbank file:
        GB_class = class_library.MakeGenbank()
        output_dict[ name ][ 'genbank_string' ] = GB_class.generate_gb_string(
            aa_seq = DB_class.aa_seq,
            fasta_name = name,
            intron_name_seq_list = kwargs[ 'intron_name_seq_list' ],
            seqlist_annotated = seqlist_annotated,
            codon2aa = DB_class.codon2aa_oneletter,
            CO_class = CO_class,
        )
        if not GB_class.check_gb(aa_seq = DB_class.aa_seq, gb_string = output_dict[ name ][ 'genbank_string' ]):
            MessageContainer.messages[name].append( '[ ERROR ] The annotation or sequence itself in the generated GenBank is not identical to the input amino acid sequence. MANUALLY VALIDATE OUTPUT!' )
        output_dict[ name ][ 'filename' ] = 'Intronserter_optDNA-{0}_{1}.gb'.format( j + 1, name.translate(remove_punctuation_map).replace(' ','')[ : 125 - 29 ] ) # remove characters that are not allowed for filenames

        DB_class.cDNA_seq_plus_i = dna_seq_fine_tuned
        output_dict[ name ][ 'cDNA_seq_plus_i' ] = DB_class.cDNA_seq_plus_i

        # Create two figures
        PlotClass = class_library.Plotting()
        output_dict[ name ][ 'fig_tmp' ] = PlotClass.plot_norm_codon_freq(dna_seq = DB_class.cDNA_seq_cleaned, aa2codon_freq_list = DB_class.aa2codon_freq_list)
        output_dict[ name ][ 'fig_tmp_introns' ] = PlotClass.plot_gene_architecture(dna_seq_list = dna_seq_list_fine_tuned, cut_sites = (kwargs[ 'cut_site_start' ], kwargs[ 'cut_site_end' ]))

        LogClass = class_library.PrepareLog()
        output_dict[ name ][ 'session_logs' ] = (
            '<ul><li>' + '</li><li>'.join([str(_) for _ in kwargs.items()]) + '</li></ul>',
            LogClass.cut_site_removal_log( log_dict = CSR_class.log_dict ),
            LogClass.intron_insertion_log( log_dict = II_class.log_dict )
            )

    kwargs[ 'output_dict' ] = output_dict
    return kwargs, MessageContainer



def get_html_strings():
    base = '''<html>
    <head>
        <script language="JavaScript">
            {functions} 
            function CopyToClipboard() {{
                const el = document.createElement('textarea');
                var text = document.getElementById("TextToCopy").innerHTML;
                el.value = text;
                document.body.appendChild(el);
                el.select();
                document.execCommand('copy');
                document.body.removeChild(el);
                alert("Copied the text: " + text);
            }}
        </script>
        <style> body {{ font-family: Arial; line-height: 1.5; }} </style>
        <style> tr:nth-child(even).tr_alternate_color {{background: #E3EBF5}} </style>
        <style> tr:nth-child(odd).tr_alternate_color {{background: #FFF}} </style>
        <style> h1, h2 {{ color: #007a00; }}  </style>
        <style> h3 {{ color: #007a00; margin-bottom: 0em; }}  </style>
    </head>'''

    html_header='''
    <body>
        <div style="word-spacing: 20px;background-color: rgba(192, 255, 33, 0.22);">
            <img src="data:image/jpeg;base64,{0}" alt="Intronserter-Logo" style="float:left;height:78px;"><br><br>
            <span style="color:white;">.</span>
            {{anchors}}
            <br><br>
        </div>
        <hr>
        <h1>Codon-optimized, cut sites removed, intron-enriched DNA sequence(s):</h1>
        {{messages}}
        <hr>'''
    html_header = html_header.format('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')

    result = '''
    <p id="{header}">
        <h2>Optimized DNA sequence of {header}:</h2>
        <table style="table-layout: fixed; width: 1024px; word-break:break-all;" >
            <tr>
                <td>>{header}</td>
            </tr>
            <tr>
                <td style="font-family:'Courier New'">{cDNA_seq_plus_i}</td>
            </tr>
        </table>
    <p>
        <table style="width: 1024px;">
            <tr>
                <td>&nbsp;&nbsp;&nbsp;&nbsp;</td>
                <td>
                    <a download="{gb_fname}" href="data:application/text;base64,{gb_base64}"><b>Download optimized DNA sequence in GenBank format</b></a>&nbsp;(should work for Firefox and Chrome, but does <u>not</u> for Internet Explorer)
                </td>
            </tr>
        </table>
    </p>
    <p>
        <table style="width: 1024px;">
            <tr>
            <td>&nbsp;&nbsp;&nbsp;&nbsp;</td>
            <td>
                <b>Please cite in your materials and methods section:</b>&nbsp;
                <i id="TextToCopy">"The sequence was optimized  based on the strategy described in Baier et al, 2018, with the tool described in Jaeger et al, 2019."</i>
                &nbsp;<a href="/chlamyintronserter?id=references">(see References)</a>
                <button onclick="CopyToClipboard()">Copy this</button>
            </td>
            </tr>
        </table>
    </p>
    <h3>Visualization of the optimization process:</h3>
    <p>
        <img src="data:image/png;base64,{fig_normfreq}">
    </p>
    <p>
        <img src="data:image/png;base64,{fig_exonintron}">
    </p>
    <button onclick="show_log_{i}()"><b>Show log of optimization process</b></button>
    <div id="logDIV{i}" style="display:none;">
        <h2>Call of Intronserter</h2>
        <p><b>Call:</b><br>{call}</p>
        <p><b>Respective content of FASTA file:</b><br>{fasta_content}</p>
        <h2>Processed Input Parameters</h2>
        <p>{params}</p>
        <h2>Cut Site Removal in codon-optimized cDNA sequence for {header}</h2>
        <p>{csr}</p>
        <h2>Intron insertion into codon-optimized and cut-site-removed cDNA sequence for {header}</h2>
        <p>{ii}</p>
        <button onclick="show_log_{i}()"><b>Hide detailed log</b></button>
    </div>
    <p><a href="#">back to top</a></p>
    <hr>
    '''

    footer = '</body></html>'

    return base, html_header, result, footer



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("aa_fasta_file", help='a file containing the FASTA AA sequence which is to be optimized', type=str)
    parser.add_argument("--output_prefix", help='prefix for the two output files (.gb and .html) (default=Intronserter_optDNA)', type=str, default='Intronserter_optDNA')
    parser.add_argument("--codon_usage_table_id", help='ID of the internally stored codon usage table (default=kazusa)', type=str, choices=['kazusa', 'hivecut'], default='kazusa')
    parser.add_argument("--custom_codon_usage_table_file", help='a file containing a codon usage table; this supersedes the parameter --codon_usage_table_id', type=str)
    parser.add_argument("--cut_sites", help='comma-separated cut sites (only 6 or 8 nt length!) to be removed from the back-translated sequence, for XbaI and XhoI e.g. TCTAGA,CTCGAG - special = custom (default=GAAGAC,GGTCTC)', type=str, default='GAAGAC,GGTCTC')
    parser.add_argument("--custom_cut_sites", help='if --cut_sites contains the entry "custom", use these additional comma-separated cut sites given as DNA sequences, e.g. TCTAGA,CTCGAG', type=str)
    parser.add_argument("--intron_seq", help='use this DNA sequence as the intron sequence (default=gtgagtcg... -> seq of rbcS2i1)', type=str, default='gtgagtcgacgagcaagcccggcggatcaggcagcgtgcttgcagatttgacttgcaacgcccgcattgtgtcgacgaaggcttttggctcctctgtcgctgtctcaagcagcatctaaccctgcgtcgccgtttccatttgcag')
    parser.add_argument("--intron_lastdifferent", help='if specified, the last intron is substituted for the seq given by --intron_lastdifferent_seq; therefore, --intron_lastdifferent_seq has be specified.', action='store_true')
    parser.add_argument("--intron_lastdifferent_seq", help='if --intron_lastdifferent is specified, use this DNA sequence as the last intron sequence', type=str)
    parser.add_argument("--supersede_intron_insert", help='if specified, no automatic determination of intron positions is performed. Instead, the positions given by --manual_intron_positions are used. Therefore, --manual_intron_positions has to be specified.', action='store_true')
    parser.add_argument("--manual_intron_positions", help='if --supersede_intron_insert is specified, use these positions instead, given as a comma-separated list, e.g. 100,450', type=str)
    parser.add_argument("--nucleotide_pair", help='nucleotide pair between which the introns are inserted (default=GG)', type=str, default='GG', choices=['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT'])
    parser.add_argument("--start", help='start exon length for automatic intron position determination (default=100)', type=int, default=100)
    parser.add_argument("--target", help='target intermediate length for automatic intron position determination (default=450)', type=int, default=450)
    parser.add_argument("--max", help='max intermediate exon length for automatic intron position determination (default=500)', type=int, default=500)
    parser.add_argument("--end", help='end exon length for automatic intron position determination (default=100)', type=int, default=100)
    parser.add_argument("--only_insert_introns", help='if specified, ONLY introns are inserted - requires a DNA sequence as input!', action="store_true")
    parser.add_argument("--cut_site_start", help='cut site to introduce at the start/5-end, e.g. None or custom or TCTAGA or CTCGAG or ... (default=None)', type=str, default='None')
    parser.add_argument("--custom_cut_site_start", help='custom cut site to introduce at the start/5-end given by this DNA sequence, e.g. TCTAGA; only active if "--cut_site_start custom" is called.', type=str)
    parser.add_argument("--cut_site_end", help='cut site to introduce at the end/3-end, e.g. None or custom or TCTAGA or CTCGAG or ... (default=None)', type=str, default='None')
    parser.add_argument("--custom_cut_site_end", help='custom cut site to introduce at the end/3-end given by this DNA sequence, e.g. TCTAGA; only active if "--cut_site_end custom" is called.', type=str)
    parser.add_argument("--linker_start", help='linker peptide to introduce at the start/5-end given by this AA sequence, e.g. GSGS (default=inactive/not set)', type=str, default='')
    parser.add_argument("--linker_end", help='linker peptide to introduce at the end/3-end given by this AA sequence, e.g. GSGS (default=inactive/not set)', type=str, default='')
    parser.add_argument("--insert_start_codon", help='if specified, a start codon is inserted at the start/5-end', action="store_true")
    parser.add_argument("--insert_stop_codon", help='if specified, a * stop codon is inserted at the end/3-end', action="store_true")
    parser.add_argument("--remove_start_codon", help='if specified, the native start codon (Met) is removed', action="store_true")
    parser.add_argument("--remove_stop_codon", help='if specified, the native * stop codon is removed', action="store_true")

    ArgsClass = parser.parse_args()
    MC_Class = MessageContainer()
    kwargs = parse_input( ArgsClass = ArgsClass, MessageContainer = MC_Class )
    kwargs, _ = process_input( kwargs = kwargs, MessageContainer = MC_Class )
    
    # display any messages
    message_list = []
    print('Info, Warning and Error messages:')
    for name, messages in MC_Class.messages.items():
        if messages:
            for message in messages:
                message_list.append( '{0}, {1}'.format(name, message) )
    if message_list:
        for message in message_list:
            print(message)
    else:
        print('None.')
    print()

    # print optimized DNA seq
    print('optimized DNA sequence(s):')
    for name, aa_seq in kwargs[ 'aa_seq_dict' ].items():
        if not aa_seq:
            continue
        print('>{0}'.format(name))
        print(kwargs[ 'output_dict' ][ name ][ 'cDNA_seq_plus_i' ])

    # save genbank file to disk
    gb_file = '{0}.gb'.format(ArgsClass.output_prefix)
    with open( gb_file, 'w' ) as fout:
        for name, aa_seq in kwargs[ 'aa_seq_dict' ].items():
            if not aa_seq:
                continue
            print( kwargs[ 'output_dict' ][ name ][ 'genbank_string' ], file = fout )

    # generate HTML file with the optimized DNA sequence and the two plots
    base, html_header, result, footer = get_html_strings()
    if not message_list:
        messages = ''
    else:
        messages= '<b><span style="color: red;">Info, Warning and Error messages:</span></b><ul>' + \
                  '<li>' + \
                  '</li><li>'.join(message_list) + \
                  '</li></ul>'
    function_template = '''            function show_log_{i}() {{
                var x = document.getElementById("logDIV{i}");
                if (x.style.display === "none") {{
                    x.style.display = "block";
                }} else {{
                    x.style.display = "none";
                }}
            }}'''

    i2fasta = {}
    with open(ArgsClass.aa_fasta_file, 'rU') as fin:
        for i, record in enumerate(SeqIO.parse(fin, "fasta")):
            i2fasta[i] = record.seq.upper()

    with open( '{0}.html'.format(ArgsClass.output_prefix), 'w' ) as fout:
        print(base.format(
                functions = os.linesep.join([function_template.format(i=i) for i in range(len(kwargs['aa_seq_dict']))])
                ),
              file=fout)
        print(html_header.format(
            anchors = '&nbsp;&nbsp;'.join(['<a href="#{header}">{header}</a>'.format(header=name) for name in kwargs['aa_seq_dict' ] ]),
            messages=messages), file=fout)
        for i, (name, aa_seq) in enumerate(kwargs[ 'aa_seq_dict' ].items()):
            if not aa_seq:
                continue
            print(
                result.format(
                    i = i,
                    header=name,
                    cDNA_seq_plus_i = kwargs[ 'output_dict' ][ name ][ 'cDNA_seq_plus_i' ],
                    gb_base64=base64.b64encode(kwargs[ 'output_dict' ][ name ][ 'genbank_string' ].encode()).decode(),
                    fig_normfreq =kwargs[ 'output_dict' ][ name ][ 'fig_tmp' ],
                    fig_exonintron =kwargs[ 'output_dict' ][ name ][ 'fig_tmp_introns' ],
                    call=sys.argv,
                    fasta_content=i2fasta[i],
                    params=kwargs[ 'output_dict' ][ name ][ 'session_logs' ][0],
                    csr=kwargs[ 'output_dict' ][ name ][ 'session_logs' ][1],
                    ii=kwargs[ 'output_dict' ][ name ][ 'session_logs' ][2],
                    gb_fname=gb_file
                ),
                file=fout
            )
        print(footer, file=fout)
