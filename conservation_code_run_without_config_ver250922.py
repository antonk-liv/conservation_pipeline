# -*- coding: utf-8 -*-
# Version 25/09/22
# Before running the code have the following files in the directory: script file, MUSCLE executable,
# proteomes.fasta, targets.fasta, sites.csv, Mapped_Uniprot_Species_Names.tsv
# Make sure you have BLAST+ installed! If code gives error, transfer all BLAST+ files into the working directory
# Enter species of target proteins using Uniprot code (see .tsv dictionary), target amino acid, likely substitution, E-value BLAST threshold in lines 20-23

###Import libraries
from Bio import SeqIO 
import os
import glob
from Bio import AlignIO
import numpy as np
import csv
import re
import pandas as pd
import shutil

###specify config settings (target species using uniprot species ID (see "Mapped_Uniprot_Species_Names.tsv"; target amino acid, substitution amino acid, Eval threshold for BLAST))
species_of_targets='HUMAN'
target_amino_acid='S'
sub='T'
eval_thres=0.00001

###create a file containing targets which are not analysed and the reason for their exclusion
excluded=open('targets_not_analysed.csv','w')
excluded.write('Protein ID'+','+'Reason for exclusion'+'\n')

### Extract target IDs from the "targets.fasta" file
target_names=[]
for targ_name in SeqIO.parse('targets.fasta', "fasta"):
    targ_NAME=targ_name.name.strip().split('|')[1]
    if targ_NAME not in target_names:
        target_names.append(targ_NAME)

### BLAST section
script_directory = os.path.dirname(os.path.realpath(__file__))
targets_input=open('targets.fasta','r') 
proteomes_input=open('proteomes.fasta','r')
#Find out the number of proteomes analysed
specs_in_proteome=[] 
for spe in SeqIO.parse('proteomes.fasta',"fasta"):
    spec=(spe.name.strip().split('|')[2]).strip().split('_')[1] #extract species name
    if spec not in specs_in_proteome:
        specs_in_proteome.append(spec)
proteome_number=int(len(specs_in_proteome)-1)
print('This analysis involves comparing', species_of_targets, 'proteins against',proteome_number, 'proteomes of other species')
print(specs_in_proteome)
#check for blast main tsv output generated from previous runs
blast_output=glob.glob(script_directory+'\*blast_result_file.tsv')
#if the output is already present, ignore BLAST step, otherwise run BLAST
if blast_output != []:
    print('BLAST output is already available from a previous run of the pipeline; BLAST step is skipped')
if blast_output == []:
    #make database with your proteomes
    command1=('makeblastdb -in %s -dbtype prot -max_file_sz 2gb -out blast_database' % 'proteomes.fasta')
    os.system(command1)
    #run your target proteins against your proteome database
    command2 = ('blastp -db blast_database -query %s -out blast_result_file.tsv -outfmt 6' % 'targets.fasta')
    os.system(command2)

###Getting top hits
tsv = pd.read_csv("blast_result_file.tsv",sep='\t',names = ["A", "B", "C", "D","E","F","G","H","I","J","K","L"])
#extract all targets
all_targets=[]
no_origin=[]
for row in range(0,len(tsv)):
    target=(str(tsv.loc[row,'A'])).strip().split('|')[1]
    species=(str(tsv.loc[row,'B'])).strip().split('|')[2] #extract species 
    species=species.strip().split('_')[1]
    if target not in all_targets and species==species_of_targets: #only extract targets for which a hit from origin species is found by BLAST
        all_targets.append(target)
#extract targets for which there are no hits from origin species        
tsv1=tsv.copy()
for roww in range(0,len(tsv1)):
    target=(str(tsv1.loc[roww,'A'])).strip().split('|')[1]
    species=(str(tsv1.loc[roww,'B'])).strip().split('|')[2]
    species=species.strip().split('_')[1]
    if target not in all_targets and target not in no_origin:
        no_origin.append(target)
        excluded.write(target+','+'No BLAST hits from target species'+'\n')
#identify targets with no matches in BLAST at all (likely due to target sequence being too short or poor annotation)
for tt in target_names:
    if tt not in all_targets  and tt not in no_origin:
        excluded.write(tt+','+'No hits at all in BLAST'+'\n')

#copy the dataframe get the top hits from each species
tsv2=tsv1.copy()
allhits_output=open("alltophits.csv","w")
allhits_output.write('Target'+','+'Species'+','+'Hit'+','+'%Identity'+','+'E-value'+','+'Max. Score'+'\n')
for t in all_targets:
    hit_output=open(t + "_hits.csv", 'w')
    hit_output.write('Target'+','+'Species'+','+'Hit'+','+'%Identity'+','+'E-value'+','+'Max. Score'+'\n')
    species_list=[] #list of all other species except humans
    hit_list=[] #list of hits from origin species
    id_list=[] #sequence simlarity of hit from origin species
    score_list=[] #BLAST score of hit from origin species
    evalues=[] #Evalue of hit from origin species
    hit_counter=0
    for rrow in range(0,len(tsv2)):
        target2=(str(tsv2.loc[rrow,'A'])).strip().split('|')[1] #extract target from origin species
        species=(str(tsv2.loc[rrow,'B'])).strip().split('|')[2] #extract species 
        species=species.strip().split('_')[1]
        top_hit=(str(tsv2.loc[rrow,'B'])).strip().split('|')[1] #extract top hit
        identity=(float(tsv2.loc[rrow,'C'])) #extract %sequence similarity
        identity=('%s' % float('%.3g' % identity)) #convert to 3 sig figs
        evalue=(str(tsv2.loc[rrow,'K'])) #extract evalue
        maxscore=(str(tsv2.loc[rrow,'L'])) #extract max score
        if t==target2 and species==species_of_targets: #extract the first hit from origin species (should be identical but sometimes it is not)
            hit_list.append(top_hit)
            id_list.append(identity)
            score_list.append(maxscore)
            evalues.append(evalue)
        if t==target2 and species==species_of_targets and t==top_hit: #what to do if the first hit from origin species was not the same as the target, but the hit which IS the same as the target is still present below (rare cases). This ensures that the top hit from origin species is always the same as the target
            hit_list[0]=top_hit
            id_list[0]=identity
            score_list[0]=maxscore
            evalues[0]=evalue
        if t==target2 and species != species_of_targets and species not in species_list and ((float(evalue)) <=eval_thres): #add top hits for the remaining species
            print(t,species,top_hit,identity,maxscore) #hit from origin species
            hit_output.write(t+','+species+','+top_hit+','+str(identity)+','+evalue+','+maxscore+'\n')
            allhits_output.write(t+','+species+','+top_hit+','+str(identity)+','+evalue+','+maxscore+'\n')
            species_list.append(species)
    if hit_list != []: #If at least 1 hit from origin species is found for the target protein (most cases) 
        print(t,species_of_targets,hit_list[0],id_list[0],evalues[0],score_list[0])
        hit_output.write(t+','+species_of_targets+','+hit_list[0]+','+str(id_list[0])+','+str(evalues[0])+','+str(score_list[0])+'\n')
        allhits_output.write(t+','+species_of_targets+','+hit_list[0]+','+str(id_list[0])+','+str(evalues[0])+','+str(score_list[0])+'\n')
    hit_output.close()
allhits_output.close()

###Getting sequences from proteome files
outputs_with_top_hits=glob.glob(script_directory+'\*_hits.csv') #searches files containing hits for each protein target
for file in outputs_with_top_hits:
    z=open(file,'r')
    hits=[] #list of all top BLAST hits for each of the species for each target
    obscure_list=[] #list of targets which are not found in the reference proteome of origin species
    true_target=[] #list containing targets which are found in the reference proteome of origin species
    
    #extract targets not in proteome, top hits, targets in proteome
    for line in z: 
        target_name=line.strip().split(',')[0]
        sp=line.strip().split(',')[1]
        hit=line.strip().split(',')[2]
        if target_name!=hit and sp==species_of_targets:
            obscure_list.append(target_name) #adds target name to the obscure list if it is not found in reference proteome (the hit from origin species is not the same as the target)
        if hit:
            hits.append(hit) #extract all hits 
        if target_name==hit and sp==species_of_targets: #extract targets which are present in the reference proteome (target=hit from origin species). This therefore eliminates targets with no hits from origin species and targets which are not in the proteome
            true_target.append(target_name)
            
    #write targets which are not in the proteome into a separate file
    if obscure_list != []:    
        for prot in obscure_list:
            excluded.write(prot+','+'Target is not the reference proteome of its species'+'\n')
         
    z.close() #close the top hits file as it is no longer needed
    hits.remove('Hit') #remove 'Hit' header from the hits list
    
    #only find hit sequences for targets which have a matching (identical) hit from the origin species
    if true_target != []: 
        for tt in true_target:
            z2=open(tt+"_hit_sequences.fasta",'w') #only create the hit sequences for those targets which are found in the reference proteome
        record_names=[]
        record_names2=[]
        for record in SeqIO.parse('proteomes.fasta', "fasta"):
            for hit in hits:
                if hit == str(((record.name).strip().split('|')[1])):
                    record_names.append(record.name)
            for element in record_names:
                if element not in record_names2: #eliminates duplicates
                    z2.write('>'+str(record.description)+'\n'+str(record.seq)+'\n')
                    record_names2.append(element)
        z2.close()
proteomes_input.close()
print('retrieving sequences for MSAs...')

###Generating multiple sequence alignments
seq_outputs=glob.glob(script_directory+'\*sequences.fasta') #searches all hit sequence files
for sequence in seq_outputs: 
    input_fasta=sequence
    input_fasta=str(input_fasta).strip().split("\\") #splits directory into sections
    input_fasta=input_fasta[len(input_fasta)-1] #identifies FASTA sequence file name in its overall directory name
    name_for_out=str(input_fasta).strip().split(".")[0] #extract file name for the alignment output
    targ_name=str(name_for_out).strip().split("_")[0] #name of target only
    all_lengths=[] #list of all sequence lengths for hits per target
    for rec in SeqIO.parse(sequence, "fasta"): #extract the lengths of all hits per target
        all_lengths.append(int(len(rec.seq)))
    max_len=int(max(all_lengths)) #identify max sequence length
    if max_len < 2000: #use default settings if max sequence length is below 2000 amino acids   
        command3 = ("muscle -in {0} -seqtype protein -out {1}.afa".format(input_fasta, name_for_out)) #default MUSCLE alignments using up to 16 iterations
        os.system(command3)
    if max_len >= 2000: #use recommended MUSCLE settings for long sequences if max sequence length is 2000 amino acids and above
        command3 = ("muscle -in {0} -seqtype protein -out {1}.afa -maxiters 2".format(input_fasta, name_for_out)) #command to run MUSCLE alignments with 2 iterations max.
        os.system(command3)
        excluded.write(targ_name+','+'Failed alignment due to one of the sequences being too long (>30000bp)'+'\n')

###Calculating conservation score of target px
alignments=glob.glob(script_directory+'\*sequences.afa') 
target_amino_acid_lowercase=target_amino_acid.lower()
conservation_file = open(target_amino_acid+'_conservation_in_'+species_of_targets+'.csv', 'w')
conservation_file.write('Protein'+','+'Peptide sequence from alignment'+','+'Position in alignment'+','+'Position in target protein'+','+'No. of species analysed'+','+'%Conservation out of total no. of proteomes'+','+'%Conservation out of total no. of proteomes inc. substitutions'+','+'%Conservation in aligned hits'+','+'%Conservation in aligned hits inc. substitutions'+','+'No. of species aligned'+','+'Species aligned (UP codes)'+','+'Species aligned (common or sci names)'+','+'No. of species aligned excl. gaps'+','+'%Conservation in aligned hits excl. gaps'+','+'%Conservation in aligned hits excl. gaps inc. substitutions'+','+'-1 site'+','+'+1 site'+','+'-1 site position in aln.'+','+'+1 site position in aln.'+','+'%Conservation of -1 site in aligned hits'+','+'%Conservation of +1 site in aligned hits'+'\n') 
aligned_orth=open('aligned_hits_per_target.csv','w') #list of all aligned hits per target
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
                #write the lists of hits matched per target into a separate CSV output
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
aligned_orth.close() #close the output with all aligned hits per target
targets_input.close()

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
#Move all top hits files into a new folder
os.mkdir('Top_BLAST_hits_per_target') #create new folder which will contain top hits for individual target protein
for top_hit_file in outputs_with_top_hits:
    shutil.move(top_hit_file, 'Top_BLAST_hits_per_target') #move each top_hits file to the folder

#Move all FASTA sequences of top hits per target to a new output folder
os.mkdir('FASTA_sequences_of_top_hits_per_target') 
for fasta_hit_file in seq_outputs:
    shutil.move(fasta_hit_file, 'FASTA_sequences_of_top_hits_per_target')

#Move all alignments for each target into a new folder
os.mkdir('Alignments_per_target')
for alignment_file in alignments:
    shutil.move(alignment_file, 'Alignments_per_target')    
    
#see if the code was run before, send the outputs to a new folder with a later number than the previous runs, making sure the code can be re-run
max_run=0
all_files=(glob.glob(script_directory+'\*'))
run_number=[] #number of runs already done
for af in all_files:
    if 'ALL_OUTPUTS_' in af:
        run_number.append(int(af[-1]))
if run_number != []: #if the code was run before, find the number of the latest run, otherwise keep as zero
    max_run=max(run_number)

#Move ALL outputs from the code into one folder
new_run_name='ALL_OUTPUTS_'+str(max_run+1)
os.mkdir(new_run_name) #create a folder which would contain all outputs
shutil.move('Top_BLAST_hits_per_target',new_run_name)
shutil.move('FASTA_sequences_of_top_hits_per_target',new_run_name)
shutil.move('Alignments_per_target',new_run_name)
shutil.move('aligned_hits_per_target.csv',new_run_name)
shutil.move('alltophits.csv',new_run_name)
#shutil.move('blast_result_file.tsv',new_run_name)
shutil.move(target_amino_acid+'_conservation_in_'+species_of_targets+'.csv',new_run_name) #change the name if target amino acid is different
#shutil.move('blast_database.phr',new_run_name)
#shutil.move('blast_database.pin',new_run_name)
#shutil.move('blast_database.psq',new_run_name)
shutil.move('targets_not_analysed.csv',new_run_name)
shutil.move('Conservation_of_target_'+target_amino_acid+'_sites.csv',new_run_name)

