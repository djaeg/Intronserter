#!/usr/bin/env python3

import os
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

    # validate the parsed fasta seq by comparison against Biopython
    with open(ArgsClass.aa_fasta_file, 'rU') as fin:
        aa_seq_dict = TMP.parse_fasta( fin.read() )
        fin.seek(0)
        for (name, seq), record in zip(aa_seq_dict.items(), SeqIO.parse(fin, "fasta")):
            #record.name == name is FALSE when the header contains any white spaces - white spaces are removed by Biopython
            if not record.seq.upper() == seq:
                MessageContainer.messages['global'].append('[ ERROR ] Parsing the FASTA AA seq of {0} was unsuccessful. VALIDATE output!'.format(name))

    # validate FASTA AA input for only valid characters
    allowed_characters = set(IUPAC.protein.letters + '*')
    assert sorted([ TMP.AA_lookup_dict[k]['1-letter'] for k in TMP.AA_lookup_dict ]) == sorted(allowed_characters)
    for name, seq in aa_seq_dict.items():
        if not set(seq).issubset(allowed_characters):
            MessageContainer.messages['global'].append('[ ERROR ] Your FASTA sequence {0} contains invalid characters. Allowed characters are {1}'.format(name, sorted(allowed_characters)))
            aa_seq_dict[name] = ''

    # check if aa- or cDNA-seq:
    allowed_characters = set('ATCG')
    for name, seq in aa_seq_dict.items():
        MessageContainer.messages[name] = []
        if seq and set(seq).issubset(allowed_characters):
            MessageContainer.messages[name].append( '[ INFO ] Your input sequence "{0}" seems to be a cDNA sequence, because it contains only the characters A, T, C and G. The sequence will be translated to its amino acid counterpart first.'.format(name) )
            if len(seq) % 3 != 0:
                MessageContainer.messages[name].append( '[ ERROR ] Your input sequence "{0}" seems to be a cDNA sequence, but its length is not a multiplier of 3! Please re-submit a valid sequence, preferably as amino acids.'.format(name) )
                aa_seq_dict[ name ] = ''
            else:
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

    kwargs[ 'output_dict' ] = None
    return kwargs



def process_input( kwargs, MessageContainer ):
    output_dict = collections.OrderedDict()
    remove_punctuation_map = dict((ord(char), None) for char in string.punctuation)

    # check if the intron sequences are free from cut sites
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
        DB_class = class_library.Tables( aa_seq.upper() )
        output_dict[ name ] = { 'cDNA_seq_plus_i' : None, 'genbank_string' : None, 'name' : '>{0}'.format(name) }

        # codon optimize:
        DB_class.import_codon_table( codon_table_string = kwargs[ 'codon_table' ] )
        DB_class.convert_codon_counts_2_freq()
        CO_class = class_library.CodonOptimizer( aa_seq = DB_class.aa_seq )
        CO_class.reverse_translate( DB_class.aa_oneletter_2_mostfreq_codon_dict )
        DB_class.cDNA_seq = CO_class.dna_seq

        # cut site removal:
        CSR_class = class_library.CutSiteRemover( dna_seq = DB_class.cDNA_seq, cut_site_list = kwargs[ 'cut_site_list' ], codon2aa = DB_class.codon2aa, aa2codon_freq_list = DB_class.aa2codon_freq_list )
        DB_class.cDNA_seq_cleaned = CSR_class.main( iter_max = 1000 )
        MessageContainer.messages[name] += CSR_class.messages
        
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
        if not CO_class.translate(coding_dna, DB_class.codon2aa_oneletter) == aa_seq or not str(Seq(coding_dna, IUPAC.unambiguous_dna).translate()) == aa_seq:
            MessageContainer.messages[name].append( '[ ERROR ] The translation of the spliced optimized DNA sequence does not match the input AA sequence. MANUALLY VALIDATE OUTPUT!' )

        # add cut sites:
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
            seqlist_annotated = seqlist_annotated
        )
        if not GB_class.check_gb(aa_seq = DB_class.aa_seq, gb_string = output_dict[ name ][ 'genbank_string' ]):
            MessageContainer.messages[name].append( '[ ERROR ] The annotation or sequence itself in the generated GenBank is not identical to the input amino acid sequence. MANUALLY VALIDATE OUTPUT!' )
        output_dict[ name ][ 'filename' ] = 'ChlamyIntronserter_optDNA-{0}_{1}.gb'.format( j + 1, name.translate(remove_punctuation_map).replace(' ','')[ : 125 - 29 ] ) # remove characters that are not allowed for filenames

        DB_class.cDNA_seq_plus_i = dna_seq_fine_tuned
        output_dict[ name ][ 'cDNA_seq_plus_i' ] = DB_class.cDNA_seq_plus_i

        # Create two figures
        PlotClass = class_library.Plotting()
        output_dict[ name ][ 'fig_tmp' ] = PlotClass.plot_norm_codon_freq(dna_seq = DB_class.cDNA_seq_cleaned, aa2codon_freq_list = DB_class.aa2codon_freq_list)
        output_dict[ name ][ 'fig_tmp_introns' ] = PlotClass.plot_gene_architecture(dna_seq_list = dna_seq_list_fine_tuned, cut_sites = (kwargs[ 'cut_site_start' ], kwargs[ 'cut_site_end' ]))

        LogClass = class_library.PrepareLog()
        output_dict[ name ][ 'session_logs' ] = (
            '<h2>Processed Input Parameters were:</h2><ul><li>' + '</li><li>'.join([str(_) for _ in kwargs.items()]) + '</li></ul>',
            '<h2>Cut Site Removal in codon-optimized cDNA sequence for {0}</h2>'.format(name) + LogClass.cut_site_removal_log( log_dict = CSR_class.log_dict ),
            '<h2>Intron insertion into codon-optimized and cut-site-removed cDNA sequence for {0}</h2>'.format(name) + LogClass.intron_insertion_log( log_dict = II_class.log_dict )
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
            <img src="data:image/jpeg;base64,{0}" alt="ChlamyIntronserter-Logo" style="float:left;width:370px;height:78px;"><br><br>
            <span style="color:white;">.</span>
            {{anchors}}
            <br><br>
        </div>
        <hr>
        <h1>Codon-optimized, cut sites removed, intron-enriched DNA sequence(s):</h1>
        {{messages}}
        <hr>'''
    html_header = html_header.format('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tHJHY+LvG2ra9YeYhXdbu0cSsM9QTCxzX1pVHQ9D0/wzothpGlWcOn6ZYwJbWtpboFjhiRQqooHQAAD8KvV7s5KUro3CiiioAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigAooooAKKKKACiiigD/9k=')

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
                    <a href="data:application/text;base64,{gb_base64}"><b>Download optimized DNA sequence in GenBank format</b></a>&nbsp;(should work for Firefox and Chrome, but does <u>not</u> for Internet Explorer)
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
        <p>{params}</p>
        <p>{csr}</p>
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
    parser.add_argument("--output_prefix", help='prefix for the two output files (.gb and .html) (default=ChlamyIntronserter_optDNA)', type=str, default='ChlamyIntronserter_optDNA')
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
                    params=kwargs[ 'output_dict' ][ name ][ 'session_logs' ][0],
                    csr=kwargs[ 'output_dict' ][ name ][ 'session_logs' ][1],
                    ii=kwargs[ 'output_dict' ][ name ][ 'session_logs' ][2],
                ),
                file=fout
            )
        print(footer, file=fout)
