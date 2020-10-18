#include "read_vcf.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<vcf_entry> load_and_filter_vcf_file( std::string input_vcf_file, int chromosome, std::string chr_string, std::string hetfilterFile ) {
        ofstream ofile_het; ofile_het.open(hetfilterFile);    // counters
        int n    = 0;  // total number of records in file
        int nsnp = 0;  // number of SNP records in file
        int nhq  = 0;  // number of SNPs for the single sample passing filters
        int nseq = 0;  // number of sequences
        int nfi_arr = 0; int nfi     = 0; int *fi     = NULL; // filter data for each call
        int ngq_arr = 0; int ngq     = 0; int *gq     = NULL; // quality data for each call
        int ndp_arr = 0; int ndp     = 0; int *dp     = NULL; // coverage data for each call
        int nad_arr = 0; int nad     = 0; int *ad     = NULL; // coverage data for each call
        int ngt_arr = 0; int ngt     = 0; int *gt     = NULL; // genotype data for each call
        std::string ref_base,var_base;
	std::list<std::string> allowed_base_list = { "A","T","C","G"};
        std::vector<vcf_entry> entry;    //kstring_t filter;
        int rec_position,rec_id,len_ref;
        int i = 0;
        htsFile * inf = bcf_open(input_vcf_file.c_str(), "r");
        bcf_hdr_t *hdr = bcf_hdr_read(inf);
        bcf_srs_t *sr = bcf_sr_init();
        bcf1_t *line;
        bcf_sr_add_reader(sr,input_vcf_file.c_str());
        while(bcf_sr_next_line(sr)) {
                line = bcf_sr_get_line(sr, 0); //cout << chromosome << "  " << line->rid << endl;
                if (line->rid == chromosome) {
                        ref_base = line->d.allele[0];
                        var_base = line->d.allele[1];
			bool ref_found = ( std::find( allowed_base_list.begin(), allowed_base_list.end(), ref_base) != allowed_base_list.end() );
			bool var_found = ( std::find( allowed_base_list.begin(), allowed_base_list.end(), var_base) != allowed_base_list.end() );
			if ( ref_found && var_found ) {
                        rec_position = line->pos;
                        rec_id = line->rid;
                        len_ref = line->rlen;  //filter = line->shared;
                        nfi = bcf_get_format_int32(hdr, line, "FI", &fi, &nfi_arr);
                        ngq = bcf_get_format_int32(hdr, line, "GQ", &gq, &ngq_arr);
                        ndp = bcf_get_format_int32(hdr, line, "DP", &dp, &ndp_arr);
                        ngt = bcf_get_format_int32(hdr, line, "GT", &gt, &ngt_arr);
                        nad = bcf_get_format_int32(hdr, line, "AD", &ad, &nad_arr);
                        int genotype = bcf_gt_allele(gt[0]);    int genotype_s = bcf_gt_allele(gt[1]);
                        //int hccbl_genotype = bcf_gt_allele(gt[2]);  int hccbl_genotype_s = bcf_gt_allele(gt[3]);
                        //int hccbl_genotype2 = bcf_gt_allele(gt[4]); int hccbl_genotype2_s = bcf_gt_allele(gt[5]);
                        if ( genotype == 0 && genotype_s == 1 ) {
                                int tot_cov = ad[0] + ad[1];
                                double fraction = std::min(ad[0]/((double)tot_cov),ad[1]/((double)tot_cov));
                                if ( ref_base.length() == 1 && var_base.length() == 1 ) {
                                if ( tot_cov > minimum_cov && tot_cov < maximum_cov ) {
                                if ( fraction > minimum_frac && fraction < (1-minimum_frac) ) {
                                        if ((DEBUG == 1) && (i > het_cutoff)) { break; }
                                        ofile_het << rec_position << "\t" << ref_base << "\t" << var_base << "\t" << fraction << "\t";
                                        ofile_het << "genotype \t" << genotype << "/" << genotype_s << "\t" << ad[0] << "\t" << ad[1] << "\t";
                                        entry.push_back(vcf_entry());
                                        entry[i].pos = rec_position;    entry[i].chromosome_id = rec_id;
                                        entry[i].ref_base = ref_base;   entry[i].var_base = var_base;
                                        i++;
                                }
                                }
                                }
                        }
			}
                }
        }
        bcf_hdr_destroy(hdr);
        bcf_close(inf);
        ofile_het.close();
        bcf_destroy(line);
        return entry;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<vcf_entry> load_vcf_file( std::string input_vcf_file, int chromosome ) {
        std::string ref_base,var_base;
        std::vector<vcf_entry> entry;
	std::list<std::string> allowed_base_list = { "A","T","C","G"};
        int rec_position,rec_id,len_ref;
        bcf_srs_t *sr = bcf_sr_init();
        bcf_sr_add_reader(sr,input_vcf_file.c_str());
        bcf1_t *line;
        int i = 0;
        while(bcf_sr_next_line(sr)) {
		line = bcf_sr_get_line(sr,0);
		if (line->rid == chromosome) {
			//if ((DEBUG == 1) && (i > het_cutoff)) { break; }
			ref_base = line->d.allele[0];
			var_base = line->d.allele[1];
                        bool ref_found = ( std::find( allowed_base_list.begin(), allowed_base_list.end(), ref_base) != allowed_base_list.end() );
                        bool var_found = ( std::find( allowed_base_list.begin(), allowed_base_list.end(), var_base) != allowed_base_list.end() ); 
			//cout << ref_base << " " << var_base << " " << ref_found << " " << var_found << endl;
                        if ( ref_found && var_found ) {
			rec_position = line->pos;	
			rec_id = line->rid;
			len_ref = line->rlen;
			entry.push_back(vcf_entry());
			entry[i].pos = rec_position;    entry[i].chromosome_id = rec_id;
			entry[i].ref_base = ref_base;   entry[i].var_base = var_base;
			entry[i].bounded = true;
			i += 1;
			}
		}
        }
      	return entry;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<vcf_entry> load_vcf_file_coverage( std::string input_vcf_file, int chromosome ) {
        std::string ref_base,var_base;
        std::vector<vcf_entry> entry;
	std::list<std::string> allowed_base_list = { "A","T","C","G"};
        int rec_position,rec_id,len_ref,nallele;
        bcf_srs_t *sr = bcf_sr_init();
        bcf_sr_add_reader(sr,input_vcf_file.c_str());
        bcf1_t *line;
        htsFile *inf = bcf_open(input_vcf_file.c_str(), "r");
        bcf_hdr_t *hdr = bcf_hdr_read(inf);
        int i = 0;
	int nformats,nsamples;
	//int nseq = 0;
	//const char **seqnames = NULL;
	int ndp_arr = 0; int ndp     = 0; int *dp     = NULL; // coverage data for each call
        int nad_arr = 0; int nad     = 0; int *ad     = NULL; // coverage data for each call
        while(bcf_sr_next_line(sr)) {
                line = bcf_sr_get_line(sr,0);
                if (line->rid == chromosome) {  //while ( bcf_read(inf, hdr, rec)>=0 ) {
                        //if ((DEBUG == 1) && (i > het_cutoff)) { break; }
			//seqnames = bcf_hdr_seqnames(hdr, &nseq);
                        ref_base = line->d.allele[0];
                        var_base = line->d.allele[1];
			bool ref_found = ( std::find( allowed_base_list.begin(), allowed_base_list.end(), ref_base) != allowed_base_list.end() );
			bool var_found = ( std::find( allowed_base_list.begin(), allowed_base_list.end(), var_base) != allowed_base_list.end() );
			if ( ref_found && var_found ) {
                        rec_position = line->pos;
                        rec_id = line->rid;
                        len_ref = line->rlen;
			nallele = line->n_allele;
			nformats = line->n_fmt;
			nsamples = line->n_sample;
                        entry.push_back(vcf_entry());
                        entry[i].pos = rec_position;    entry[i].chromosome_id = rec_id;
                        entry[i].ref_base = ref_base;   entry[i].var_base = var_base;
                        entry[i].bounded = true;
			ndp = bcf_get_format_int32(hdr, line, "DP", &dp, &ndp_arr);
			nad = bcf_get_format_int32(hdr, line, "AD", &ad, &nad_arr);
			uint8_t *ptr = (uint8_t*)line->indiv.s;
			bcf_fmt_t fmt[6];
			cout << rec_position << " " << nallele << " " << rec_id << " " << nformats << " " << nsamples << endl;
    			for (int k=0; k<nallele; k++) { cout << k << " " << line->d.allele[k] << " " << ad[k] << " " << endl; }
			//for (int k=0; k<nformats; k++) { ptr = bcf_unpack_fmt_core1(ptr, nsamples, &fmt[k]); }
			//for (int k=0; k<nformats; k++) { ptr = bcf_unpack_fmt_core1(ptr, nsamples, fmt[k]); }
			for (int k=0; k<bcf_hdr_nsamples(hdr); k++){
			}
                        i += 1; //}
			}
                }
        }
        bcf_hdr_destroy(hdr);
        bcf_close(inf);
        bcf_destroy(line);
        return entry;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<vcf_entry> load_vcf_file_total( std::string input_vcf_file ) {
        std::string ref_base,var_base;
        std::vector<vcf_entry> entry;
	std::list<std::string> allowed_base_list = { "A","T","C","G"};
        int rec_position,rec_id,len_ref;
        bcf_srs_t *sr = bcf_sr_init();
        bcf_sr_add_reader(sr,input_vcf_file.c_str());
        bcf1_t *line;
        int i = 0;
        while(bcf_sr_next_line(sr)) {
                line = bcf_sr_get_line(sr,0);
                ref_base = line->d.allele[0];
                var_base = line->d.allele[1];
		bool ref_found = ( std::find( allowed_base_list.begin(), allowed_base_list.end(), ref_base) != allowed_base_list.end() );
		bool var_found = ( std::find( allowed_base_list.begin(), allowed_base_list.end(), var_base) != allowed_base_list.end() );
		if ( ref_found && var_found ) {
                rec_position = line->pos;
                rec_id = line->rid;
                len_ref = line->rlen;
                entry.push_back(vcf_entry());
                entry[i].pos = rec_position;    entry[i].chromosome_id = rec_id;
                entry[i].ref_base = ref_base;   entry[i].var_base = var_base;
                entry[i].bounded = true;
                i += 1;
		}
        }
        return entry;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
std::map<std::string,int> load_vcf_file_header( std::string input_vcf_file, std::map<std::string,int>& chr_str_map ) {
        htsFile *test_bcf = NULL;
        bcf_hdr_t *test_header = NULL;
        test_bcf = bcf_open(input_vcf_file.c_str(), "r");
        if(test_bcf == NULL) { throw std::runtime_error("Unable to open file."); }
        test_header = bcf_hdr_read(test_bcf);
        if(test_header == NULL) { throw std::runtime_error("Unable to read header."); }
	int nseq = 0;
	const char **seqnames = NULL;
	seqnames = bcf_hdr_seqnames(test_header, &nseq);
	for (int i = 0; i < nseq; i++) { std::string id_string = bcf_hdr_id2name(test_header,i); chr_str_map[id_string] = i; }
        bcf_hdr_destroy(test_header);
        bcf_close(test_bcf);
	return chr_str_map;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
void subset_het_sites( vcf_vector& vvec, int start_bound, int end_bound ) {
	for (int i=0; i < vvec.size(); i++) {
		if ((vvec[i].pos <= start_bound) || (vvec[i].pos >= end_bound)) {
			vvec[i].bounded = false;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////



        //bcf1_t *test_record = bcf_init();
        //bcf_destroy(test_record); 
	//cout << id_string << "\t" << i << endl;

        //std::cout << "chromosome\tposition\tnum_alleles" << std::endl;

		//if (i > 30) { break; }
        //while(bcf_read(test_bcf, test_header, test_record) == 0) {
        //    std::string id_string = bcf_hdr_id2name(test_header, test_record->rid);
        //    std::cout << id_string << "\t" << test_record->rid << "\t" << test_record->pos << "\t" << test_record->n_allele << "\t" << std::endl;
        //}


////////////////////////////////////////////////////////////////////////////////////////////////////////
//std::vector<vcf_entry> load_and_filter_vcf_file( std::string input_vcf_file, int chromosome, std::string chr_string, std::string hetfilterFile ) {
//        ofstream ofile_het; ofile_het.open(hetfilterFile);    // counters
//        int n    = 0;  // total number of records in file
//        int nsnp = 0;  // number of SNP records in file
//        int nhq  = 0;  // number of SNPs for the single sample passing filters
//        int nseq = 0;  // number of sequences
//        int nfi_arr = 0; int nfi     = 0; int *fi     = NULL; // filter data for each call
//        int ngq_arr = 0; int ngq     = 0; int *gq     = NULL; // quality data for each call
//        int ndp_arr = 0; int ndp     = 0; int *dp     = NULL; // coverage data for each call
//        int nad_arr = 0; int nad     = 0; int *ad     = NULL; // coverage data for each call
//        int ngt_arr = 0; int ngt     = 0; int *gt     = NULL; // genotype data for each call
//        std::string ref_base,var_base;
//        std::vector<vcf_entry> entry;    //kstring_t filter;
//        int rec_position,rec_id,len_ref;
//        int i = 0;
//        htsFile * inf = bcf_open(input_vcf_file.c_str(), "r");
//        bcf_hdr_t *hdr = bcf_hdr_read(inf);
//        bcf_srs_t *sr = bcf_sr_init();
//        bcf1_t *line;
//        bcf_sr_add_reader(sr,input_vcf_file.c_str());
//        while(bcf_sr_next_line(sr)) {
//                line = bcf_sr_get_line(sr, 0); //cout << chromosome << "  " << line->rid << endl;
//                if (line->rid == chromosome) {
//                        ref_base = line->d.allele[0];
//                        var_base = line->d.allele[1];
//                        rec_position = line->pos;
//                        rec_id = line->rid;
//                        len_ref = line->rlen;  //filter = line->shared;
//                        nfi = bcf_get_format_int32(hdr, line, "FI", &fi, &nfi_arr);
//                        ngq = bcf_get_format_int32(hdr, line, "GQ", &gq, &ngq_arr);
//                        ndp = bcf_get_format_int32(hdr, line, "DP", &dp, &ndp_arr);
//                        ngt = bcf_get_format_int32(hdr, line, "GT", &gt, &ngt_arr);
//                        nad = bcf_get_format_int32(hdr, line, "AD", &ad, &nad_arr);
//                        int hcc_genotype = bcf_gt_allele(gt[0]);    int hcc_genotype_s = bcf_gt_allele(gt[1]);
//                        int hccbl_genotype = bcf_gt_allele(gt[2]);  int hccbl_genotype_s = bcf_gt_allele(gt[3]);
//                        int hccbl_genotype2 = bcf_gt_allele(gt[4]); int hccbl_genotype2_s = bcf_gt_allele(gt[5]);
//                        if ( hccbl_genotype2 == 0 && hccbl_genotype2_s == 1 ) {
//                                int tot_cov = ad[4] + ad[5];
//                                double fraction = std::min(ad[4]/((double)tot_cov),ad[5]/((double)tot_cov));
//                                if ( ref_base.length() == 1 && var_base.length() == 1 ) {
//                                if ( tot_cov > minimum_cov && tot_cov < maximum_cov ) {
//                                if ( fraction > minimum_frac && fraction < (1-minimum_frac) ) {
//                                        if ((DEBUG == 1) && (i > het_cutoff)) { break; }
//                                        ofile_het << rec_position << "\t" << ref_base << "\t" << var_base << "\t" << fraction << "\t";
//                                        ofile_het << "HCC1954 \t" << hcc_genotype << "/" << hcc_genotype_s << "\t" << ad[0] << "\t" << ad[1] << "\t";
//                                        ofile_het << "HCC1954 BL \t" << hccbl_genotype << "/" << hccbl_genotype_s << "\t" << ad[2] << "\t" << ad[3] << "\t";
//                                        ofile_het << "HCC1954BL \t" << hccbl_genotype2 << "/" << hccbl_genotype2_s << "\t" << ad[4] << "\t" << ad[5] << endl;
//                                        entry.push_back(vcf_entry());
//                                        entry[i].pos = rec_position;    entry[i].chromosome_id = rec_id;
//                                        entry[i].ref_base = ref_base;   entry[i].var_base = var_base;
//                                        i++;
//                                }
//                                }
//                                }
//                        }
//                }
//        }
//        bcf_hdr_destroy(hdr);
//        bcf_close(inf);
//        ofile_het.close();
//        bcf_destroy(line);
//        return entry;
//}
