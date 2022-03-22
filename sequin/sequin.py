#!/usr/bin/env python3

import pydna.all
import Bio
import dna_features_viewer
import pathlib
import sequin.config as config
import os

class record:
    def __init__(self, oligos=dict(), sequences=dict()):
        '''Create a new record. Oligos and sequences dicts can be provided already.'''
        self.oligos = oligos
        self.sequences = sequences
        self.folder = pathlib.Path(config.folder)
        self.prefix = config.prefix

    def __repr__(self):
        return f'Clown record with {len(self.oligos)} oligos, {len(self.sequences)} sequences.'

    def fragmentize(self, name):
        '''Figure out if a parameter is a fragment name or a raw fragment variable, and returns the fragment.'''
        sequence = self.sequences[name] if type(name) == str else name
        return sequence
     
    def primerize(self, oligo):
        pr = pydna.primer.Primer

        if type(oligo) == str:
            if oligo in self.oligos.keys():
                return pr(self.oligos[oligo]) # Called by name
            else:
                return pr(oligo) # Direct sequence input
        else:
            return pr(oligo) # Just in case

    ### Plotting
    def to_graphic_record(self, sequence, topology):
        '''Make a graphic record from a sequence.'''
        translator = dna_features_viewer.BiopythonTranslator()
        translator.default_feature_color = config.default_color
        
        graphic_record = translator.translate_record(sequence, topology)
        return graphic_record

    def show_map(self, name, zoom=None, topology=None, width=config.default_plot_width):
        '''Display a map of a DNA fragment in Dseqrecord format.'''
        sequence = self.sequences[name] if type(name) == str else name
        if not topology:
            if zoom or not sequence.circular:
                topology = 'linear'
            else:
                topology = 'circular'
        
        graphic_record = self.to_graphic_record(sequence, topology)
        
        if zoom:
            graphic_record = graphic_record.crop(zoom)
        
        if len(graphic_record.features) >= config.max_features:
            print (f'Too many features ({len(sequence.features)}) to be plotted. Try cropping or increasing the maximum.')
            return
        
        ax,_ = graphic_record.plot(figure_width=width)
    
    def show_sequence(self, name, location, 
                      width=config.default_plot_width,
                      highlight=None, translate=None, strand=+1):
        '''Zoom on a part of the sequence with the "location" coordinates.'''
        sequence = self.sequences[name] if type(name) == str else name

        if location[1]-location[0] > config.max_sequence_length:
            print (f'Sequence too long to be displayed.')
            return
    
        if strand == -1:
            sequence = sequence.reverse_complement()
            l = len(sequence)
            location = (l-location[1],l-location[0])
            if highlight:
                highlight = (l-highlight[1],l-highlight[0])
            if translate:
                translate = (l-translate[1],l-translate[0])

        translator = dna_features_viewer.BiopythonTranslator()
        translator.default_feature_color = config.default_color
        graphic_record = translator.translate_record(sequence,'linear').crop(location)
        ax,_ = graphic_record.plot(figure_width=width, plot_sequence=True)

        if highlight:
            ax.fill_between([highlight[0]-0.5, highlight[1]-0.5], -0.5, -.9,
                            facecolor="orange", alpha=0.5)
        if translate:
            # Trim the interval to be a multiple of 3
            excess = translate[1]-translate[0]%3 
            if excess > 0:
                translate = (translate[0], translate[1]-excess)
            graphic_record.plot_translation(ax, translate,
                    fontdict={'weight': 'bold'}, long_form_translation=False)
        
    #### I/O
    def find_gb(self, filename, folder=None, extension='.gb', name=None):
        '''Find a gb file that matches the name (with .gb extension), and import it as a sequence.'''
        if not folder:
            folder = self.folder

        if not name:
            name = filename

        for root, dirs, files in os.walk(folder):
            if filename+extension in files:
                path = os.path.join(root, filename+extension)
                print(f'Adding {path} as {name}.')
                return self.add_gb(path, name=name)
        print(f'No matching file found with extension {extension}')

    def add_gb(self, filename, name=None, circular=True):
        '''Import a gb file to the sequence dict.'''
        seq = next(Bio.SeqIO.parse(filename, 'gb'))
        dseq = pydna.all.Dseqrecord(seq, circular=circular)
        if name:
            self.sequences[name] = dseq
        else:
            return dseq

    def write_gb(self, name, filename=None, folder=None):
        '''Write a sequence to a gb file.'''
        if not filename:
            filename = name
        fullname = str(pathlib.Path(folder) / filename) + '.gb'

        if not folder:
            folder = self.folder

        fragment = self.fragmentize(name)
        print(type(fragment))
        if type(fragment) in [pydna.contig.Contig, pydna.dseqrecord.Dseqrecord]:
            if fragment.locus == 'name':
                fragment.locus = filename if filename else name # Set the locus variable to the fragment's name
            fragment.write(fullname, f='gb')
        else:
            print ('Input should be Dseqrecord or String')
    
    def add_oligos(self, filename, sep='\t', header=True, name_col=0, seq_col=1):
        '''Parse oligos from a tsv format'''
        with open(filename,'r') as data:
            if header:
                next(data) # Skip colnames line
            count,existing = 0,0
            for line in data:
                if line != '\n':
                    ls = line.split(sep)
                    if ls[0] in self.oligos.keys():
                        existing += 1
                    else:
                        count += 1
                        self.oligos[ls[name_col].strip()] = ls[seq_col].strip()
            print(f'Found new {count} oligos ({existing} existing).')
    
    def subset(self, target, keys):
        '''Returns the subset of a dictionary with only the specified keys.'''
        return {key:value for key,value in target.items() if key in keys}

    #### Features
    def clean_features(self, name, verbose = True,
                       junk=['name','source'],
                       junk_types=['primer','primer_bind']):
        '''Remove features whose labels is in the junk list.'''
        sequence = self.sequences[name] if type(name)==str else name
        feats = sequence.features
        cleaned = []
        for i in feats:
            keep = True
            if not 'label' in i.qualifiers.keys():
                if verbose: print(f'Removing feature "{i}" that has no label.')
                keep = False
            elif i.qualifiers['label'][0] in junk:
                if verbose: print(f'Removing junk feature "{i.qualifiers["label"][0]}".')
                keep = False
            elif i.type in junk_types:
                if verbose: print(f'Removing primer {i.qualifiers["label"][0]}')
                keep = False

            if keep:
                cleaned.append(i)
        if type(name)==str:
            self.sequences[name].features = cleaned
        else:
            sequence.features = cleaned
            return sequence
            
    def new_feature(self, sequence, location, strand, name, feature_type, force=False):
        '''Add a custom feature to a sequence.'''
        feature_location = Bio.SeqFeature.FeatureLocation(location[0], location[1], strand)
        new_feature = Bio.SeqFeature.SeqFeature(feature_location, type=feature_type)
        new_feature.qualifiers['label'] = [name]

        existing_locations = [i.location for i in self.sequences[sequence].features]
        if (not force) and (feature_location in existing_locations):
            print('An identical feature is already present at this location. Skipping.')
            return
        else:
            self.sequences[sequence].features.append(new_feature)

    def list_features(self, name):
        '''List features of a sequence.'''
        sequence = self.sequences[name] if type(name) == str else name
        for i in sequence.features:
            if 'label' in i.qualifiers.keys():
                print (f'{i.qualifiers["label"][0]}: {i.location} ({i.type})')
            else:
                print (i.type, i.qualifiers)

    def find_features(self, name, targets, margins=(0,0), plot=True):
        '''Find coordinates of a region englobing all matching features. margins is a tuple giving how many base pairs to include before and after. Return a tuple of coordinates.'''
        sequence = self.sequences[name] if type(name) == str else name
        targets = [targets] if type(targets) == str else targets
        
        features = [i for i in sequence.features if 'label' in i.qualifiers.keys()]
        match = []
        for i in features:
            for j in targets:
                label = i.qualifiers['label'][0]
                if j in label:
                    print(label + ': ' + str(i.location))
                    match.append(i)

        if not match:
            print('No feature found!')
            return

        start = min([i.location.start for i in match])-margins[0]
        end = max([i.location.end for i in match])+margins[1]
        
        if plot:
            self.show_map(sequence, zoom=(start,end))
        return (start, end)
    
    def find_motif(self, target, motif, strand):
        '''Find all occurences of a given DNA motif in a DNA strand.'''
        sequence = self.sequences[target] if type(target) == str else target
        # Reverse the motif if necessary
        if strand==-1: motif = Bio.Seq.Seq(motif).reverse_complement()

        # Format the strings correctly
        motif = str(motif).upper()
        sequence = str(sequence.seq).upper()

        # Do the search
        return Bio.SeqUtils.nt_search(sequence, motif)[1:]


    def find_orfs(self, target, strand, min_length=300, start_codon='ATG', greedy=True):
        '''Detect ORFs in the strand. Returns the first and last nucleotides' coordinates.
        Greedy mode will only return the longest ORF associated with each stop codon.'''
        sequence = self.sequences[target] if type(target) == str else target
        
        # Check strand
        if strand not in [-1,1]:
            raise Exception("Strand must be either -1 and +1")

        if strand==-1: sequence = sequence.reverse_complement()

        # Find start codons
        start_codons = self.find_motif(sequence, start_codon, +1)

        orfs = []
        # Find ORFs
        for i in start_codons:
            endpoint = (len(sequence)-i-1)%3 + 1
            protein = sequence[i:-endpoint].translate()
            stop = i+(protein.seq.find('*')+1)*3
            if stop-i >= min_length:
                orfs.append((i,stop))
        
        # Filter for the longest ORFs if greedy is enabled
        if greedy:
            stops = set([i[1] for i in orfs])
            longest = {stop:max([i[1]-i[0] for i in orfs if i[1]==stop]) for stop in stops}
            orfs = [i for i in orfs if i[1]-i[0]==longest[i[1]]]

        # Reverse the orfs if looking at the -1 strand
        if strand == -1:
            l = len(sequence)
            orfs = [(l-i[1], l-i[0]) for i in orfs]

        return orfs

    def find_primers(self, target, oligos=None,
                     plot=True, zoom=None, 
                     lim=config.min_primer_length):
        '''Find primers that bind on a given target sequence. If oligos are not specified, use the record's oligo dict. Otherwise, oligos can be provided as a dict, a list/tuple, or a single oligo.'''
        sequence = self.sequences[target] if type(target) == str else target
        
        # If oligos are not specified, used the stored dict
        oligos = self.oligos if not oligos else oligos
            
        # Convert to dict
        if type(oligos) != dict:
            if type(oligos) in [list,tuple]:
                oligos = {('oligo_'+str(n)):i for n,i in enumerate(oligos)}
            else:
                oligos = {'oligo':oligos}

        # Make a list of oligos in pydna-friendly format
        oligos_pydna = []
        for i in oligos.keys():
            o = pydna.primer.Primer(oligos[i])
            o.name = i
            oligos_pydna.append(o)

        # Find binding sites
        anneal = pydna.all.Anneal(template=sequence, 
                                  primers=oligos_pydna,
                                  limit=lim)

        # Plot the map with the primers on top
        if plot:
            ft_size = len(sequence)*0.1
            gr = self.to_graphic_record(sequence, 'linear')

            gf = dna_features_viewer.GraphicFeature
            for i in anneal.forward_primers:
                gr.features.append(gf(start=i.position, end=i.position-len(i.footprint),
                                     strand=-1, color='#90bc21', label=i.name))
            for i in anneal.reverse_primers:
                gr.features.append(gf(start=i.position+len(i.footprint), end=i.position,
                                     strand=+1, color='#c65423', label=i.name))
            if zoom: gr = gr.crop(zoom)
            gr.plot()

        return anneal.report()


    #### Basic operations
    def pcr(self, template, F, R, name=None,
            verbose=config.verbose, lim=config.min_primer_length):
        '''Simulate a polymerase chain reaction.'''
        sequence = self.fragmentize(template)
       
        pr = pydna.primer.Primer
        oF = self.primerize(F)
        oR = self.primerize(R)
        oF.name = F
        oR.name = R
        amplicon = pydna.all.pcr(oF, oR, sequence, limit=lim)
        
        tmF = self.tm(amplicon.forward_primer.footprint)
        tmR = self.tm(amplicon.reverse_primer.footprint)
        
        if name and type(template)==str:
            print(f'{name}\t{template}\t{F}\t{R}\t{len(amplicon)}/{tmF}/{tmR}')
        else:
            print(f'{template} + {F} + {R} = {len(amplicon)} bp, {tmF}-{tmR}°C.')
        
        # Remove primer features
        self.clean_features(amplicon, junk=[], junk_types=['primer_bind'], verbose=False) 

        if name:
            self.sequences[name] = amplicon
        else:
            return amplicon

    def pair(self, o1, o2):
        '''Anneal together a pair of oligos, either called by their name or directly input as a sequence.'''
        if type(o1) == str:
            o1_seq = self.oligos[o1] if o1 in self.oligos.keys() else o1
        else:
            o1_seq = str(o1)
        if type(o2) == str:
            o2_seq = self.oligos[o2] if o2 in self.oligos.keys() else o2
        else:
            o2_seq = str(o2)

        return pydna.all.Dseqrecord(pydna.all.Dseq(o1_seq, o2_seq))


    def digest(self, target, enzymes, fragment=0, name=None):
        rb = Bio.Restriction.RestrictionBatch(enzymes)
        product = self.sequences[target].cut(rb)

        print(f'Digestion yields {len(product)} fragments:')
        for n,i in enumerate(product):
            extracted = '(extracted)' if n==fragment else ''
            print (f'{n}: {len(i)} bp {extracted}')
        
        if name:
            self.sequences[name] = product[fragment]
        else:
            return product[fragment]

    def point_mutation(self, template, mutation, name=None):
        '''Introduces a point mutation in a sequence. Input is a string: <before><position><after>.'''
        sequence = self.sequences[template] if type(template) == str else template
        
        if type(mutation) != str:
            print ('`mutation` should be a string: <before><position><after>.')

        # Parse the mutation string
        position = int(mutation[1:-1])-1
        before =  mutation[0].upper()
        after = mutation[-1].upper()
        
        # Make a mutated version of the sequence
        mutated = Bio.Seq.Seq(str(sequence.seq)).tomutable()

        # Introduce the mutation
        original = mutated[position].upper()
        if original != before:
            print (f'The nucleotide at {position} is {original}, not {before}.')
        mutated[position] = after

        # Reconstruct the seqrecord with features
        sequence_mutated = pydna.all.Dseqrecord(record=mutated.toseq(), 
                                                circular=sequence.circular) 
        sequence_mutated.features = sequence.features
        if name:
            self.sequences[name] = sequence_mutated
        else:
            return sequence_mutated

    #### Assemblies
    def check_homology(self, s1, s2, minimum=10, maximum=100, show_linker=True):
        '''Check that two sequences have homology regions suitable for Gibson assembly.'''
        for i in range(minimum,maximum):
            hom1 = str(s1.seq[-i-1:-1]).lower()
            hom2 = str(s2.seq[0:i]).lower()
            if hom1==hom2:
                if show_linker:
                    print(hom1)
                return i+1
        return 0

    def assemble(self, fragments, name=None,
               lim=config.min_primer_length, verbose=config.verbose, skip_check=False):
        '''Simulate isothermal homology-based assemblies like Gibson's.'''
        fragments_seq = []
        for i in fragments:
            fragment = self.sequences[i] if type(i)==str else i
            fragments_seq.append(fragment)

        if not skip_check:
            # Check that the fragments can be assembled together
            junctions = '//- '
            for i in range(len(fragments_seq)-1):
                hom = self.check_homology(fragments_seq[i], fragments_seq[i+1])
                junctions += f'{len(fragments_seq[i])} bp - [{hom} bp] -'
            hom = self.check_homology(fragments_seq[-1], fragments_seq[0])
            junctions += f'{len(fragments_seq[-1])} bp - [{hom} bp] -//'
            print (junctions)

        asm = pydna.all.Assembly(fragments_seq, limit=lim)
        assembled = asm.assemble_circular()
        if len(assembled) > 0:
            product = assembled[0]
        else:
            print ('Assembly failed!')
            return (asm)
            
        print(f'{name}: {len(product)} bp circular product.')
        if verbose:
            print(product.figure())

        if name:
            self.sequences[name] = product
        else:
            return product
        
    def ligate(self, fragments, name=None):
        '''Simulate ligation. Fragments can be either sequences or names.'''
        if type(fragments) != list:
            print ('Fragments must be supplied as a list of either names or sequences')
            return

        fragments_seq = [self.fragmentize(i) for i in fragments]
        
        # Add reversed versions of all the fragments
        fragments_seq_rev = [i.reverse_complement() for i in fragments_seq]
        fragments_u = fragments_seq + fragments_seq_rev

        # List all 5' ends
        five = [i.seq.five_prime_end()[1] for i in fragments_u]
        # List all 3' ends
        three = [pydna.utils.rc(i.seq.three_prime_end()[1]) for i in fragments_u]

        # Put the fragments in the right order
        graph = [0]
        while True:
            last = three[graph[-1]]
            if last in five: # Move on to the next fragment
                next_fragment = five.index(last)
                if next_fragment == 0:
                    circular = True
                    break
                else:
                    graph.append(next_fragment)
                    five[next_fragment] = None
            else:
                circular = False
                break
        print(f'Fragments assembled in this order: {graph}. ({"circular" if circular else "linear"})')

        product = fragments_u[0] # Reconstruct the assembled sequence
        for i in graph[1:]: 
            product += fragments_u[i]
        if circular: # Circularize if needed
            product = product.looped()
        if name:
            self.sequences[name] = product
        else:
            return product

    #### Design
    def tm(self, primer):
        '''Estimates the melting point of a primer.'''
        return round(pydna.tm.tm_default(primer),2)
    
    def make_primers(self, target, location, names=None, tm=58, lim=18,
            plot=True, plot_annealing=False, print_tm=False):
        '''Creates a pair of primers to amplify the region `location` from target.'''
        sequence = self.sequences[target] if type(target) == str else target
        
        start = location[0]
        end = location[1]
        template = sequence[start:end]

        amplicon = pydna.all.primer_design(template, target_tm=tm, limit=lim)
        f = amplicon.forward_primer
        r = amplicon.reverse_primer
        
        if plot:
            self.show_map(amplicon)
        if plot_annealing:
            self.show_sequence(sequence, (max(0,start-20),start+len(r)+20), 
                               highlight=(start,start+len(f)))
            self.show_sequence(sequence, (end-len(r)-20, min(len(sequence),end+20)),
                               highlight=(end-len(r),end))
        if print_tm:
            print(f'F: {f.seq} (Tm={self.tm(f.seq)})')
            print(f'R: {r.seq} (Tm={self.tm(r.seq)})')

        if names:
            print(f'{names[0]}\t{f.seq}\n{names[1]}\t{r.seq}')
            return {names[0]:f.seq, names[1]:r.seq}
        else:
            return f.seq,r.seq

    def make_junction(self, left, right, linker=None, plot=True, homology=40):
        '''Design tailed oligos for homology-based assembly. Input is 1) the reverse primer of the first fragment, 2) the forward primer of the second fragment, 3) optionally a linker to insert between the two fragments.'''
        if not linker:
            linker = ''

        top = Bio.Seq.reverse_complement(left).lower() + linker.lower() + right.upper()
        bottom = Bio.Seq.reverse_complement(right).lower() + Bio.Seq.reverse_complement(linker).lower() + left.upper()

        # Trim the oligos if they are longer than necessary, and equilibrate lengths
        trim = max(0,round((len(top)-homology)/2))
        top_oligo = top[trim:]
        bottom_oligo = bottom[trim:]

        # Show the junction
        print ('Junction:')
        print ("5' " + " "*trim + top_oligo + " 3'")
        print ("3' " + bottom_oligo[::-1] + " "*trim + " 5'")

        return (top_oligo, bottom_oligo)

    def latest(self, prefix=None):
        '''Looks at already existing oligos in the form `prefix+integer`, and finds what the next one should be.'''
        if not prefix:
            prefix = config.prefix

        indices = [int(i[len(prefix):]) for i in self.oligos.keys() if i[0:len(prefix)]==prefix]
        return prefix + str(max(indices)+1)

    def new_oligos(self, oligos, prefix=None):
        '''Determines if each oligo from the list is new. If it is, finds a name for it in the form `prefix+integer`.'''
        if not prefix:
            prefix = self.prefix

        oligos = oligos if type(oligos) in [list,tuple] else [oligos]
        new = dict()

        existing_names = list(self.oligos.keys())
        existing_seqs = [str(j).upper() for j in list(self.oligos.values())]

        for i in oligos:
            if str(i).upper() in existing_seqs:
                index = existing_seqs.index(str(i).upper())
                name = existing_names[index]
                print (f'There is already an oligo ({name}) with sequence {i}.')
            else:
                name = self.latest(prefix)
                self.oligos[name] = i
            # Make a dictionary of the oligos with their names
            new[name] = i
        
        print ('\nName\tSequence')
        for name,seq in sorted(new.items()):
            print (f'{name}\t{seq}') # Copy-pastable format
        return new



