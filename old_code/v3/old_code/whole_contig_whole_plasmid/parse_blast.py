all_hits_dict['canu'] = {}
canu_hits_dict = all_hits_dict['canu']

for file in canu_blast_results:

    path = Path(file)
    suffix = path.suffixes[0]
    assembly_id = path.name.split(suffix)[0].replace("_canu_blast", "")

    if assembly_id not in canu_hits_dict:
        canu_hits_dict[assembly_id] = {}

    current_barcode = canu_hits_dict[assembly_id]

    res_handle = open(file)
    blast_records = NCBIXML.parse(res_handle)

    for record in blast_records:
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                if (hsp.expect == 0 or hsp.expect < 1*10^-100):
                    query_length = record.query_length
                    contig_id = record.query.split(" ")[0]
                    hit_id = alignment.title.split(" ")[1]
                    alignment_length = alignment.length
                    bit_score = hsp.score
                    evalue = hsp.expect
                    if contig_id not in current_barcode:
                        current_barcode[contig_id] = [{'ID' : hit_id, 'alignment_length' : alignment_length, 'query_length' : query_length, 'bit_score' : bit_score, 'e-value' : evalue}]
                    else:
                        current_barcode[contig_id].append({'ID' : hit_id, 'alignment_length' : alignment_length, 'query_length' : query_length, 'bit_score' : bit_score, 'e-value' : evalue})
                    #print("****Alignment****")
                    #print("Barcode: ", assembly_id)
                    #print("contig ID: ", contig_id)
                    #print("sequence:", hit_id)
                    #print("length:", alignment.length)
                    #print("bit score:",hsp.score)
                    #print("e value:", hsp.expect)

pprint.pprint(all_hits_dict['canu'])
#out_file = open(output_json+'canu_hits.json', "w")
#json.dump(all_hits_dict['canu'],out_file, indent=6)