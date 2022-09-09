# -*- coding: utf-8 -*-
#The code allows to process alignments (make sure they are in the same folder as this code) and identify conservation of target amino acid at every position in sequence
#Make sure file "Mapped_Uniprot_Species_Names.tsv" is in the directory, as well as the alignments and corresponding files with complete FASTA sequences of proteins involved in alignments
#Make sure that the file containing the positions of target sites and corresponding proteins (as UniProt Accession IDs) is also included in the working directory 
#Make sure that configurations files "configurations_aln_file.ini" and "link_to_config_aln.py" are in the working directory
#Make sure that the following modules are installed:
from Bio import SeqIO 
import os
import glob
from Bio import AlignIO
import numpy as np
import re
import csv
import shutil

###specify config settings in configurations_aln_file.ini file (target species using uniprot species ID (see "Mapped_Uniprot_Species_Names.tsv"; target amino acid, substitution amino acid, number of proteomes of interest)
import link_to_config_aln
config = link_to_config_aln.read_config()
species_of_targets = config['Conservation_Code_Parameters_Aln_Section']['species_of_targets']
target_amino_acid = config['Conservation_Code_Parameters_Aln_Section']['target_amino_acid']
sub = config['Conservation_Code_Parameters_Aln_Section']['sub']
proteome_number = float(config['Conservation_Code_Parameters_Aln_Section']['proteome_number'])

###directory
script_directory = os.path.dirname(os.path.realpath(__file__))

#Grab all alignments in .afa format and calculate conservation of target amino acid
alignments=glob.glob(script_directory+'\*sequences.afa') 

###create a file containing targets which are not analysed and the reason for their exclusion
excluded=open('targets_not_analysed.csv','w')
excluded.write('Protein ID'+','+'Reason for exclusion'+'\n')

target_amino_acid_lowercase=target_amino_acid.lower()
conservation_file = open(target_amino_acid+'_conservation_in_'+species_of_targets+'.csv', 'w')
conservation_file.write('Protein'+','+'Peptide sequence from alignment'+','+'Position in alignment'+','+'Position in target protein'+','+'No. of species analysed'+','+'%Conservation out of total no. of proteomes'+','+'%Conservation out of total no. of proteomes inc. substitutions'+','+'%Conservation in aligned orthologues'+','+'%Conservation in aligned orthologues inc. substitutions'+','+'No. of species aligned'+','+'Species aligned (UP codes)'+','+'Species aligned (common or sci names)'+','+'No. of species aligned excl. gaps'+','+'%Conservation in aligned orthologues excl. gaps'+','+'%Conservation in aligned orthologues excl. gaps inc. substitutions'+','+'-1 site'+','+'+1 site'+','+'-1 site position in aln.'+','+'+1 site position in aln.'+','+'%Conservation of -1 site in aligned orthologues'+','+'%Conservation of +1 site in aligned orthologues'+'\n') 
aligned_orth=open('aligned_orthologues_per_target.csv','w') #list of all aligned orthologues per target
aligned_orth.write('Target protein'+','+'No. of aligned species'+','+'Species aligned (UP codes)'+','+'Species aligned (common or sci names)'+'\n')
for file in alignments:
    #process alignment
    if os.path.getsize(file) != 0: #if the alignment is not empty
        f=open(file,'r')
        prot_NAME=str(file) #directory of your alignment file
        NAME=prot_NAME.strip().split('\\') #directory split into section by \
        len_NAME=(len(NAME))-1 #index of the alignment
        NAME=prot_NAME.strip().split('\\')[int(len_NAME)] #extract the name of target protein from the directory by its index and species
        name_final=NAME.strip().split('_')[0] #extract target protein name only
        
        alignment=AlignIO.read(f, "fasta")
        species_names=[] #list of all species names present in the alignment
        record_names=[]
        target_protein_names=[]
        human_seq='' #target sequence (blank to begin with but will be extracted later)
        for record in alignment:
            if species_of_targets in record.id:
                target_name_from_aln=((record.id).strip().split('|')[1]) #extracts name of target protein from record
                target_protein_names.append(target_name_from_aln)
                human_seq=str(record.seq)
                    
            if record.id not in record_names:
                hit_name=((record.id).strip().split('|')[2])
                species_name=((hit_name).strip().split('_')[1])
                species_names.append(str(species_name)) #extracts species names
                record_names.append(record.id)
        
        if len(species_names) > 1: #only process alignments with more than 1 species (target matched with at least 1 other sequence)
            if target_amino_acid not in human_seq: #eliminate targets with no target amino acid in their sequence
                excluded.write(target_name_from_aln+','+'No target amino acid in sequence'+'\n')
                f.close() #close alignment
            else:
                #map target names in alignments (given as Uniprot mnemonic codenames) to common or scientific species names
                names_to_map=[]
                for item in species_names:
                    if item != species_of_targets:
                        names_to_map.append(item)
                mnemonics_from_UP=[] #all mnemonic codenames from Uniprot dictionary
                common_names_from_UP=[] #all common names from UP
                sci_names_from_UP=[] #all scientific names from UP
                uniprot_vocab=open("Mapped_Uniprot_Species_Names.tsv",'r')
                for line in uniprot_vocab: #loop through the UP dictionary and extract all data
                    mnemonic=line.strip().split('\t')[0]
                    mnemonics_from_UP.append(mnemonic)
                    common_name=line.strip().split('\t')[1]
                    common_names_from_UP.append(common_name)
                    sci_name=line.strip().split('\t')[2]
                    sci_names_from_UP.append(sci_name)
                uniprot_vocab.close()
                UP_array=np.stack((mnemonics_from_UP,common_names_from_UP,sci_names_from_UP),axis=1) #array of all name mappings
                mapped_names=[] #final list of mapped names
                for codename in names_to_map:
                    if codename not in UP_array: #if target codename is not found in UP dictionary then map it to itself
                        mapped_names.append(codename)
                    else:
                        for term in UP_array:
                            if codename == term[0]: #map codename to common species name
                                mapped_name=term[1]
                                if mapped_name=='NA': #map codename to sci name if no common name is found
                                    mapped_name=term[2]
                                mapped_names.append(mapped_name)
                #convert the lists of menmonic names and mapped names to single strings (names separated by ;) to allow to add them into the output
                mapped_names_as_string=''
                for mn in mapped_names:
                    mapped_names_as_string+=mn
                    mapped_names_as_string+=';' #convert all mapped names to a single string (can't use comma due to CSV file being used as output)
                mnemonics_as_string=''
                for m in names_to_map:
                    mnemonics_as_string+=m
                    mnemonics_as_string+=';'
                #write the lists of orthologues matched per target into a separate CSV output
                aligned_orth.write(target_name_from_aln+','+str(len(mapped_names))+','+mnemonics_as_string+','+mapped_names_as_string+'\n')
                
                #process the alignment to retrieve target sequence, positions of target amino acid and conservation scores
                human_pos_in_alignment=int(species_names.index(species_of_targets)) #specifies row position of target sequence in the alignment
                number_of_seqs=len(alignment)
                alignment_length=len(alignment[0].seq)
                aa_positions=[] #lists positions in the alignment where target amino acid from a target protein is found (counting from 0)
                aa_positions2=[] #lists positions in the alignment (counting from 1)
                aa_positions_minus=[] #lists -1 positions from the target amino acid
                aa_positions_plus=[] #lists +1 positions from the target amino acid
                amino_acids_at_pos=[] #all amino acids across sequences at position where target amino acid from a target sequence is found
                amino_acids_at_minus_pos=[] #all amino acids at -1 position next to the target amino acid from target species
                amino_acids_at_plus_pos=[] #all amino acids at +1 position next to the target amino acid from target species
                minus_pos_in_aln=[] #actual positions of -1 sites in aln
                plus_pos_in_aln=[] #actual positions of +1 sites in aln
                conservation_scores=[] #a list of conservation scores for each target amino acid within target protein
                conservation_scores2=[] #a list of conservation scores for target protein including interchangeable amino acids
                conservation_scores3=[] #a list of conservation scores for targets out of the number of proteomes analysed
                conservation_scores3_1=[] #a list of conservation scores for targets out of the number of proteomes analysed incl. mutations
                conservation_scores_nogaps=[] #a list of conservation scores for target amino acid(excludes a species from total count if gap is found)
                conservation_scores2_nogaps=[] #a list of conservation scores for target including interchangeable amino acids (excluding gaps)
                conservation_scores_minus=[] #a list of conservation scores of -1 sites next to a target site
                conservation_scores_plus=[] #a list of conservation scores of +1 sites next to a target site
                # get the positions of where target amino acid from target sequence is in the alignment
                for j in range(0, alignment_length): #loop through each column in the alignment
                    chars_at_pos = [] #amino acid residues at a column in the alignment
                    for i in range(0,number_of_seqs): #loop through each row per column and add an amino acid in each row to the list
                        chars_at_pos.append(alignment[i].seq[j]) 
                    if target_amino_acid in chars_at_pos[human_pos_in_alignment]: #find column positions containing a target amino acid in the target sequence
                            alignment_pos=(int(j))
                            aa_positions.append(int(j)) #finds alignment positions where target amino acid is present and adds to list of positions
                            aa_positions2.append(int(j)+1)
                for column in aa_positions:
                    amino_acids_at_pos.append(alignment[:, column]) #finds amino acids across all sequences at a certain position if target amino acid is found at this position
                #get the -1 and +1 positions around target amino acid within target sequence
                for poz in aa_positions:
                    #-1 sites:
                    number=1 #numbers to be added to -1/+1 shifts if there is a gap where a proximal site is meant to be
                    poz_minus=poz-number #original proximal position: target amino acid minus 1
                    if poz_minus!=int(-1): #only add the -1 position if it exists (it's possible that alignment begins with target aa and no -1 site can therefore be present)
                        while (alignment[:, poz_minus][human_pos_in_alignment]) == '-' and poz_minus >0: #keep shifting left if there is a gap in alignment where -1 site is meant to be until the actual site is found. Also making sure that the shifting does not go below position 0
                            number+=1
                            poz_minus=poz-number
                        if (alignment[:, poz_minus][human_pos_in_alignment]) == '-': #if there are only gaps at -1 position to the target
                            aa_positions_minus.append('Not found')
                            minus_pos_in_aln.append('Not found')
                        else:
                            aa_positions_minus.append(poz_minus)
                            minus_pos_in_aln.append(poz_minus+1)   
                    if poz_minus == int(-1): #if there is no -1 position and the alignment begins with target amino acid write not found
                        aa_positions_minus.append('Not found') 
                        minus_pos_in_aln.append('Not found')
                    #+1 sites
                    number2=1
                    poz_plus=poz+number2                
                    if poz_plus!=alignment_length: #only add +1 position if it doesn't exceed alignment length (can happen if target sequence aa is the last amino acid)
                        while (alignment[:, poz_plus][human_pos_in_alignment]) == '-' and poz_plus < alignment_length-1: #keep shifting left if there is a gap in alignment where -1 site is meant to be until the actual site is found. If last position is reached and there is a gap, still add the gap
                            number2+=1
                            poz_plus=poz+number2
                        if (alignment[:, poz_plus][human_pos_in_alignment]) == '-': #if there are only gaps on the +1 side of the target, write 'Not found' for proximal site info 
                            aa_positions_plus.append('Not found')
                            plus_pos_in_aln.append('Not found')
                        else:
                            aa_positions_plus.append(poz_plus)
                            plus_pos_in_aln.append(poz_plus+1)
                    if poz_plus==alignment_length:
                        aa_positions_plus.append('Not found')
                        plus_pos_in_aln.append('Not found')     
                #get the actual amino acids at -1 and +1 positions in the alignment the same way it was done for the target
                for column2 in aa_positions_minus:
                    if column2 != 'Not found':    
                        amino_acids_at_minus_pos.append(alignment[:, column2])
                    else: 
                        amino_acids_at_minus_pos.append('Not found')
                for column3 in aa_positions_plus:
                    if column3 != 'Not found':    
                        amino_acids_at_plus_pos.append(alignment[:, column3])
                    else: 
                        amino_acids_at_plus_pos.append('Not found')
                #Extracting 15aa peptide sequences containing every target residue as well as proximal +1 and -1 sites around the target. All positional scenarios are included
                peptide_seq=[] #lists all peptide sequences containing target amino acid across all species in the alignment
                species_total_nogaps_list=[]
                species_total_list=[]
                seq_file=open(name_final+'_hit_sequences.fasta')
                records = list(SeqIO.parse(seq_file, "fasta"))
                final_15aa_sequences=[]
                left_sites=[]
                right_sites=[]
                for record in records:
                    if name_final in record.id:
                        positions_in_protein_code=[]
                        positions_in_protein=[] #list of actual target amino acid positions in the protein sequence
                        record_seq=str(record.seq)
                        pattern=target_amino_acid
                        regex = re.compile(pattern)
                        for match in regex.finditer(record_seq):
                            posi=match.start()
                            lower_aa=record_seq[posi].lower()
                            record_seq=record_seq[:posi] + record_seq[posi + 1:] 
                            record_seq=record_seq[:posi]+lower_aa+record_seq[posi:]
                            pos_in_protein=posi+1 #actual position of the target site in the protein sequence
                            positions_in_protein.append(pos_in_protein)
                            positions_in_protein_code.append(posi)
                        for prot_pos in positions_in_protein_code:
                            if int(prot_pos)==6:
                                final_sequence=record_seq[prot_pos-6:prot_pos+9]
                                left=record_seq[prot_pos-1]
                                left_sites.append(left)
                                right=record_seq[prot_pos+1]
                                right_sites.append(right)
                                if len(final_sequence)==15:
                                    final_15aa_sequences.append(final_sequence)
                            if int(prot_pos)==5:
                                final_sequence=record_seq[prot_pos-5:prot_pos+10]
                                left=record_seq[prot_pos-1]
                                left_sites.append(left)
                                right=record_seq[prot_pos+1]
                                right_sites.append(right)
                                if len(final_sequence)==15:
                                    final_15aa_sequences.append(final_sequence)
                            if int(prot_pos)==4:
                                final_sequence=record_seq[prot_pos-4:prot_pos+11]
                                left=record_seq[prot_pos-1]
                                left_sites.append(left)
                                right=record_seq[prot_pos+1]
                                right_sites.append(right)
                                if len(final_sequence)==15:
                                    final_15aa_sequences.append(final_sequence)
                            if int(prot_pos)==3:
                                final_sequence=record_seq[prot_pos-3:prot_pos+12]
                                left=record_seq[prot_pos-1]
                                left_sites.append(left)
                                right=record_seq[prot_pos+1]
                                right_sites.append(right)
                                if len(final_sequence)==15:
                                    final_15aa_sequences.append(final_sequence)
                            if int(prot_pos)==2:
                                final_sequence=record_seq[prot_pos-2:prot_pos+13]
                                left=record_seq[prot_pos-1]
                                left_sites.append(left)
                                right=record_seq[prot_pos+1]
                                right_sites.append(right)
                                if len(final_sequence)==15:
                                    final_15aa_sequences.append(final_sequence)
                            if int(prot_pos)==1:
                                final_sequence=record_seq[prot_pos-1:prot_pos+14]
                                left=record_seq[prot_pos-1]
                                left_sites.append(left)
                                right=record_seq[prot_pos+1]
                                right_sites.append(right)
                                if len(final_sequence)==15:
                                    final_15aa_sequences.append(final_sequence)
                            if int(prot_pos)==0:
                                final_sequence=record_seq[prot_pos:prot_pos+15]
                                right=record_seq[prot_pos+1]
                                right_sites.append(right)
                                left='Not found'
                                left_sites.append(left)
                                if len(final_sequence)==15:
                                    final_15aa_sequences.append(final_sequence)
                            if int(prot_pos)==len(record_seq):
                                final_sequence=record_seq[prot_pos-15:prot_pos]
                                left=record_seq[prot_pos-1]
                                left_sites.append(left)
                                right='Not found'
                                right_sites.append(right)
                                if len(final_sequence)==15:
                                    final_15aa_sequences.append(final_sequence)
                            if int(prot_pos)==(len(record_seq)-1):
                                final_sequence=record_seq[prot_pos-14:prot_pos+1]
                                left=record_seq[prot_pos-1]
                                left_sites.append(left)
                                right='Not found'
                                right_sites.append(right)
                                if len(final_sequence)==15:
                                    final_15aa_sequences.append(final_sequence)
                            if int(prot_pos)==(len(record_seq)-2):
                                final_sequence=record_seq[prot_pos-13:prot_pos+2]
                                left=record_seq[prot_pos-1]
                                left_sites.append(left)
                                right=record_seq[prot_pos+1]
                                right_sites.append(right)
                                if len(final_sequence)==15:
                                    final_15aa_sequences.append(final_sequence)
                            if int(prot_pos)==(len(record_seq)-3):
                                final_sequence=record_seq[prot_pos-12:prot_pos+3]
                                left=record_seq[prot_pos-1]
                                left_sites.append(left)
                                right=record_seq[prot_pos+1]
                                right_sites.append(right)
                                if len(final_sequence)==15:
                                    final_15aa_sequences.append(final_sequence)
                            if int(prot_pos)==(len(record_seq)-4):
                                final_sequence=record_seq[prot_pos-11:prot_pos+4]
                                left=record_seq[prot_pos-1]
                                left_sites.append(left)
                                right=record_seq[prot_pos+1]
                                right_sites.append(right)
                                if len(final_sequence)==15:
                                    final_15aa_sequences.append(final_sequence)
                            if int(prot_pos)==(len(record_seq)-5):
                                final_sequence=record_seq[prot_pos-10:prot_pos+5]
                                left=record_seq[prot_pos-1]
                                left_sites.append(left)
                                right=record_seq[prot_pos+1]
                                right_sites.append(right)
                                if len(final_sequence)==15:
                                    final_15aa_sequences.append(final_sequence)
                            if int(prot_pos)==(len(record_seq)-6):
                                final_sequence=record_seq[prot_pos-9:prot_pos+6]
                                left=record_seq[prot_pos-1]
                                left_sites.append(left)
                                right=record_seq[prot_pos+1]
                                right_sites.append(right)
                                if len(final_sequence)==15:
                                    final_15aa_sequences.append(final_sequence)
                            if int(prot_pos)==(len(record_seq)-7):
                                final_sequence=record_seq[prot_pos-8:prot_pos+7]
                                left=record_seq[prot_pos-1]
                                left_sites.append(left)
                                right=record_seq[prot_pos+1]
                                right_sites.append(right)
                                if len(final_sequence)==15:
                                    final_15aa_sequences.append(final_sequence)
                            if int(prot_pos) > 6 and int(prot_pos) < (len(record_seq)-7): #all positions at which you can have +7 and -7 amino acids around the target amino acid
                                final_sequence=record_seq[prot_pos-7:prot_pos+8]
                                left=record_seq[prot_pos-1]
                                left_sites.append(left)
                                right=record_seq[prot_pos+1]
                                right_sites.append(right)
                                if len(final_sequence)==15:
                                    final_15aa_sequences.append(final_sequence)
                #find % conservation at every position where a target amino acid from target species is found and set the score to 0 if no target amino acid is found in any other species
                for amino_acid in amino_acids_at_pos:
                    aa_count=amino_acid.count((target_amino_acid))-1 #count of target amino acid excluding origin species
                    aa_count2=aa_count+amino_acid.count(sub) #count of target + count of interchangeable amino acid in relation to target aa
                    species_total=(int(len(amino_acid)))-1 #total number of sequences (denominator including gaps) minus origin species
                    species_total_list.append(species_total)
                    gap_count=amino_acid.count('-') #count of gaps at position
                    species_total_nogaps=species_total-gap_count #total number of sequences (denominator excluding gaps)
                    species_total_nogaps_list.append(species_total_nogaps)
                    print(target_name_from_aln,species_total,species_total_nogaps,amino_acid, aa_count)
                    if aa_count==0 and aa_count2==0: #if there are no target amino acids or substitutions found at the position of the target amino acid from target species make the cons. score 0
                        conservation_score=0
                        conservation_score2=0
                        conservation_score3=0
                        conservation_score3_1=0
                        conservation_score_nogaps=0
                        conservation_score2_nogaps=0
                        conservation_scores.append(conservation_score)
                        conservation_scores2.append(conservation_score)
                        conservation_scores3.append(conservation_score3)
                        conservation_scores3_1.append(conservation_score3_1)
                        conservation_scores_nogaps.append(conservation_score_nogaps)
                        conservation_scores2_nogaps.append(conservation_score2_nogaps)
                    if aa_count==0 and aa_count2!=0: #if there are no target amino acids but only substitutions, for example, for target S there are only T residues in alignment
                        conservation_score=0
                        conservation_score2=('%s' % float('%.3g' % (int(aa_count2)/int(species_total)*100))) #conservation score of target aa including interchangeable amino acids in the numerator
                        conservation_score3=0
                        conservation_score3_1=('%s' % float('%.3g' % (int(aa_count2)/(proteome_number)*100)))
                        conservation_score_nogaps=0
                        conservation_score2_nogaps=('%s' % float('%.3g' % (int(aa_count2)/int(species_total_nogaps)*100))) 
                        conservation_scores.append(conservation_score)
                        conservation_scores2.append(conservation_score)
                        conservation_scores3.append(conservation_score3)
                        conservation_scores3_1.append(conservation_score3_1)
                        conservation_scores_nogaps.append(conservation_score_nogaps)
                        conservation_scores2_nogaps.append(conservation_score2_nogaps)
                    if aa_count!=0: #calculate conservation score if there is at least 1 target amino acid found at the target position in the alignment in other species
                        conservation_score=('%s' % float('%.3g' % (int(aa_count)/int(species_total)*100))) #generates conservation score to 3 sig. figs.
                        conservation_score2=('%s' % float('%.3g' % (int(aa_count2)/int(species_total)*100))) #conservation score at target position including interchangeable amino acids in the numerator
                        conservation_score3=('%s' % float('%.3g' % (int(aa_count)/(proteome_number)*100)))
                        conservation_score3_1=('%s' % float('%.3g' % (int(aa_count2)/(proteome_number)*100)))
                        conservation_score_nogaps=('%s' % float('%.3g' % (int(aa_count)/int(species_total_nogaps)*100)))
                        conservation_score2_nogaps=('%s' % float('%.3g' % (int(aa_count2)/int(species_total_nogaps)*100)))
                        conservation_scores.append(conservation_score)
                        conservation_scores2.append(conservation_score2)
                        conservation_scores3.append(conservation_score3)
                        conservation_scores3_1.append(conservation_score3_1)
                        conservation_scores_nogaps.append(conservation_score_nogaps)
                        conservation_scores2_nogaps.append(conservation_score2_nogaps)
                #get coservation scores of -1 and +1 sites around target. Consider scenarios where the site is a gap in alignment and where there are 0 matches in the alignment (site from origin species present only)
                #-1 conservation
                for site in amino_acids_at_minus_pos:
                    if site=='Not found': #notes if no -1 site is found
                        conservation_score_minus='Not found'
                        conservation_scores_minus.append(conservation_score_minus)
                    else:
                        minus_count=site.count(site[human_pos_in_alignment])-1 #Count of -1 amino acids at a -1 position in alignment excluding the target amino acid (-1)
                        if minus_count==0: #make conservation score 0 if all other sites in the alignment at the -1 position are gaps or the target site from target sequence is unique in the alignment
                            conservation_score_minus=0
                            conservation_scores_minus.append(conservation_score_minus)
                        if minus_count!=0 and site!='Not found':
                            conservation_score_minus=('%s' % float('%.3g' % (int(minus_count)/int(species_total)*100)))
                            conservation_scores_minus.append(conservation_score_minus)
                #+1 conservation
                for site2 in amino_acids_at_plus_pos:
                    if site2=='Not found': #if no amino acid is found at the +1 position (the sequence ends with target amino acid)
                        conservation_score_plus='Not found'
                        conservation_scores_plus.append(conservation_score_plus)
                    else: #if proximal site is found and it's not a gap, get its conservation score
                        plus_count=site2.count(site2[human_pos_in_alignment])-1 #The count of +1 amino acids at a +1 target position in alignment excluding the target amino acid (-1)
                        if plus_count==0:
                            conservation_score_plus=0
                            conservation_scores_plus.append(conservation_score_plus)
                        if plus_count!=0 and site2!='Not found':
                            conservation_score_plus=('%s' % float('%.3g' % (int(plus_count)/int(species_total)*100)))
                            conservation_scores_plus.append(conservation_score_plus)
                #Combine ALL the information above into one array and print in the output
                seq_pos_score_array=np.stack((final_15aa_sequences,aa_positions2,positions_in_protein,conservation_scores3,conservation_scores3_1,conservation_scores,conservation_scores2,species_total_list,conservation_scores_nogaps,conservation_scores2_nogaps,species_total_nogaps_list,left_sites,right_sites,minus_pos_in_aln,plus_pos_in_aln,conservation_scores_minus,conservation_scores_plus),axis=1) #3D array of target peptide sequence containing target amino acid, its position in the alignment and %conservation    
                for value in seq_pos_score_array:
                    conservation_file.write(target_name_from_aln+','+value[0]+','+value[1]+','+value[2]+','+str(proteome_number)+','+value[3]+','+value[4]+','+value[5]+','+value[6]+','+value[7]+','+mnemonics_as_string+','+mapped_names_as_string+','+value[10]+','+value[8]+','+value[9]+','+value[11]+','+value[12]+','+value[13]+','+value[14]+','+value[15]+','+value[16]+'\n')
                f.close() #close alignment
                seq_file.close() #close file with FASTA sequences
        else: #if only 1 hit is found (target itself), add to the specific file
            excluded.write(target_name_from_aln+','+'Target has no significant hits in BLAST'+'\n')
            f.close() #close alignment

#close outputs
excluded.close()
conservation_file.close() #close output containing peptide sequences with target px
aligned_orth.close() #close the output with all aligned orthologues per target

###Cross-reference conservation data with sites.csv file containing positions of target sites for which conservation needs to be analysed
file1=open(target_amino_acid+'_conservation_in_'+species_of_targets+'.csv','r') #file from which data will be extracted (target file)
file2=open('sites.csv','r') #Template file into which target data will be added
file1_reader=csv.reader(file1,delimiter=',')
file2_reader=csv.reader(file2,delimiter=',')
header_f1=next(file1_reader) #ignore header
header_f2=next(file2_reader)
out_file  = open('Conservation_of_target_'+target_amino_acid+'_sites.csv', "w") #output
#include the header into output from template file and then add necessarry headings from target file to the header in template
for h in header_f2:
    out_file.write(h+',')
for hh in header_f1:
    out_file.write(hh+',')
out_file.write('\n')
#Extract all the needed data from target file to be added to template file
all_target_data=[] #each line in target file
all_target_prots=[] #each target in target file
for line in file1_reader:
    if line[0] not in all_target_prots:
        all_target_prots.append(line[0]) #extract prot name from target data
    all_target_data.append(line)
#extract all data from the template file
all_template_data=[] #each line in template file
for line2 in file2_reader:
    all_template_data.append(line2)
###cross-reference two files
found=[] #all cross-referenced lines
for x in all_template_data: #for each LINE in template data
    #identify columns in template file which need to be cross-referenced with target file
    temp_poz=x[1]
    temp_prot=x[0]
    if ';' in temp_prot:
        temp_prot=temp_prot.strip().split(';')[0]
    #cross-reference those columns with target file if template protein is found in target data
    if temp_prot in all_target_prots:
        for y in all_target_data: #for each LINE in target data
            targ_prot=y[0]
            targ_poz=y[3]
            if temp_prot==targ_prot and temp_poz==targ_poz:
                for col in x:
                    out_file.write(col+',')
                #write in all columns from target file separated by comma
                for yyy in y:
                    out_file.write(str(yyy)+',')
                out_file.write('\n')
                found.append(x)
    #Any line in the template file not found in the target will be written into the output and dashes will be added to the columns where target data should have been
    if x not in found:
        x += len(header_f1) * ['-']
        for coll in x:
            out_file.write(coll+',')
        out_file.write('\n')               
file1.close()
file2.close()
out_file.close()

###Move all individual outputs per target to new corresponding folders:

#Move all FASTA sequences of top hits per target to a new output folder
seq_outputs=glob.glob(script_directory+'\*sequences.fasta')
os.mkdir('FASTA_sequences_of_top_hits_per_target') 
for fasta_hit_file in seq_outputs:
    shutil.move(fasta_hit_file, 'FASTA_sequences_of_top_hits_per_target')

#Move all alignments for each target into a new folder
os.mkdir('Alignments_per_target')
for alignment_file in alignments:
    shutil.move(alignment_file, 'Alignments_per_target')    
    
#Move ALL outputs from the code into one folder
os.mkdir('ALL_OUTPUTS') #create a folder which would contain all outputs
shutil.move('FASTA_sequences_of_top_hits_per_target','ALL_OUTPUTS')
shutil.move('Alignments_per_target','ALL_OUTPUTS')
shutil.move(target_amino_acid+'_conservation_in_'+species_of_targets+'.csv','ALL_OUTPUTS') #change the name if target amino acid is different
shutil.move('targets_not_analysed.csv','ALL_OUTPUTS')
shutil.move('Conservation_of_target_'+target_amino_acid+'_sites.csv','ALL_OUTPUTS')


