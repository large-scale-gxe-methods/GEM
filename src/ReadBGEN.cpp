
#include "declars.h"
#include "ReadBGEN.h"





void Bgen::processBgenHeaderBlock(char genofile[300]) {


	// Ensure genotype file contains BGEN extension
	string genopath(genofile);
	if (genopath.substr(genopath.length() - 5, 5) != ".bgen") {
		cout << genopath << " is not a .bgen file. Currently only supporting genotype files in BGEN format. \n";
		exit(1);
	}


	/***************************************************************
	  Read General information of bgen data.
	  Conduct Sample IDMatching process if necessary.
	***************************************************************/
	fin = fopen(genofile, "rb");
	if (fin == 0) {
		cerr << "BGEN file could not be opened." << endl;
		exit(1);
	}



	cout << "General information of BGEN file. \n";


	// First four bytes (offset)
	fread(&offset, 4, 1, fin); cout << "offset: " << offset << '\n';
	

	// The header block
	uint L_H;      fread(&L_H, 4, 1, fin);   cout << "L_H: \n";
	fread(&Mbgen, 4, 1, fin); cout << "BGEN snpBlocks (Mbgen): " << Mbgen << '\n'; assert(Mbgen != 0);
	fread(&Nbgen, 4, 1, fin); cout << "BGEN samples (Nbgen): " << Nbgen << '\n';
	char magic[5]; fread(magic, 1, 4, fin); magic[4] = '\0';  //cout << "magic bytes: " << string(magic) << endl;
	fseek(fin, L_H - 20, SEEK_CUR);                           //cout << "skipping L_H-20 = " << L_H-20 << " bytes (free data area)" << endl;
	uint flags;    fread(&flags, 4, 1, fin);                  //cout << "flags: " << flags << endl;


	// The header block - flag definitions
	CompressedSNPBlocks = flags & 3; cout << "CompressedSNPBlocks: " << CompressedSNPBlocks << '\n';
	assert(CompressedSNPBlocks == 1); // REQUIRE CompressedSNPBlocks==1
	Layout = (flags >> 2) & 0xf; cout << "Layout: " << Layout << '\n';
	assert(Layout == 1 || Layout == 2); // REQUIRE Layout==1 or Layout==2
	SampleIdentifiers = flags >> 31; cout << "SampleIdentifiers: " << SampleIdentifiers << '\n';


}




std::vector<long int> getPositionOfBgenVariant(FILE* fin, uint offset, uint Mbgen, uint Nbgen, uint CompressedSNPBlocks, uint Layout, vector<int> Mbgen_begin) {

	int t = 0;
	uint maxLA = 65536;
	uint maxLB = 65536;
	char* snpID   = new char[maxLA + 1];
	char* rsID    = new char[maxLA + 1];
	char* chrStr  = new char[maxLA + 1];
	char* allele1 = new char[maxLA + 1];
	char* allele0 = new char[maxLB + 1];

	vector <uchar> zBuf;
	vector <long int> variant_pos(Mbgen_begin.size());



	fseek(fin, offset + 4, SEEK_SET);



	for (int snploop = 0; snploop < Mbgen; snploop++) {

		 if (snploop == Mbgen_begin[t]) {
			 variant_pos[t] = ftell(fin);
			 t++;
		 }
		 
		

		 /**** Variant Data Block ********/

		 // Number of individuals. Only present when Layout == 1
		 if (Layout == 1) {
			 uint Nrow; fread(&Nrow, 4, 1, fin); // cout << "Nrow: " << Nrow << " " << std::flush;  
			 if (Nrow != Nbgen) {
				 cerr << "ERROR: Nrow = " << Nrow << " does not match Nbgen = " << Nbgen << '\n';
				 exit(1);
			 }
		 }

		 // The length of the variant identifier
		 ushort LS; fread(&LS, 2, 1, fin);  // cout << "LS: " << LS << " " << std::flush;
		 if (LS > maxLA) {
			 maxLA = 2 * LS;
			 delete[] snpID;
			 char* snpID = new char[maxLA + 1];
		 }
		


		 // The variant identifier
		 fread(snpID, 1, LS, fin); snpID[LS] = '\0'; // cout << "snpID: " << string(snpID) << " " << std::flush;
		

		 // The length of the rsid
		 ushort LR; fread(&LR, 2, 1, fin); // cout << "LR: " << LR << " " << std::flush;
		 if (LR > maxLA) {
			 maxLA = 2 * LR;
			 delete[] rsID;
			 char* rsID = new char[maxLA + 1];
		 }

		 // The rsid
		 fread(rsID, 1, LR, fin); rsID[LR] = '\0'; // cout << "rsID: " << string(rsID) << " " << std::flush;


		 // The length of the chromosome
		 ushort LC; fread(&LC, 2, 1, fin); // cout << "LC: " << LC << " " << std::flush;

		 // The chromosome
		 fread(chrStr, 1, LC, fin); chrStr[LC] = '\0';

		 // The variant position
		 uint physpos; fread(&physpos, 4, 1, fin); // cout << "physpos: " << physpos << " " << std::flush;


		 // The number of alleles if Layout == 2. If Layout == 1, this value is assumed to be 2
		 if (Layout == 2) {
			 ushort LKnum; fread(&LKnum, 2, 1, fin); // this is for Layout = 2, Lnum = 2 is Layout = 1

			 if (LKnum != 2) {
				 cerr << "ERROR: Non-bi-allelic variant found: " << LKnum << " alleles\n";
				 exit(1);
			 }
		 }


		 // Length of the first allele
		 uint LA; fread(&LA, 4, 1, fin); // cout << "LA: " << LA << " " << std::flush;
		 if (LA > maxLA) {
			 maxLA = 2 * LA;
			 delete[] allele1;
			 char* allele1 = new char[maxLA + 1];
		 }
		 // The first allele
		 fread(allele1, 1, LA, fin); allele1[LA] = '\0';


		 // The length of the second allele
		 uint LB; fread(&LB, 4, 1, fin); // cout << "LB: " << LB << " " << std::flush;
		 if (LB > maxLB) {
			 maxLB = 2 * LB;
			 delete[] allele0;
			 char* allele0 = new char[maxLB + 1];
		 }
		 // The second allele
		 fread(allele0, 1, LB, fin); allele0[LB] = '\0';


		 /* Reading genotype for BGEN v1.2/3 */
		 uint zLen;  fread(&zLen, 4, 1, fin); // cout << "zLen: " << zLen << endl;
		 uint DLen;
		 if (CompressedSNPBlocks == 0) {
			 fseek(fin, zLen, SEEK_CUR);

		 } else {
			 fseek(fin, 4 + zLen - 4, SEEK_CUR);
		 }

	}


	return variant_pos;
}














void BgenParallelGWAS(int begin, int end, long int byte, char genobgen[300], int thread_num,  Bgen test) {


	boost::thread_specific_ptr<BgenThreadSpecificVariables> tss;
	//boost::thread_specific_ptr<vector <double>> AF_vec;
	if (!tss.get()) {
		tss.reset(new BgenThreadSpecificVariables);
	}

	std::string output = "gem_bin_" + std::to_string(thread_num) + ".tmp";
	std::ofstream results(output, std::ofstream::binary);

	int samSize = test.samSize;
	int numSelCol = test.numSelCol;
	uint CompressedSNPBlocks = test.CompressedSNPBlocks;
	int robust = test.robust;
	int Sq = test.Sq;
	double* resid2 = test.resid;
	vector<long int> include_idx = test.include_idx;
	(*tss).stream_snps = test.stream_snps;
	(*tss).ZGSvec.resize(samSize * (1 + Sq) * (*tss).stream_snps);
	(*tss).ZGSR2vec.resize(samSize * (1 + Sq) * (*tss).stream_snps);
	(*tss).WZGSvec.resize(samSize * (1 + Sq) * (*tss).stream_snps);
	(*tss).WZGS = &(*tss).WZGSvec[0];
	(*tss).Sq1 = Sq + 1;
	(*tss).geno_snpid.resize((*tss).stream_snps);
	(*tss).fin3 = fopen(genobgen, "rb");


	fseek((*tss).fin3, byte, SEEK_SET);



	(*tss).snploop = begin;
	while ((*tss).snploop <= end) {

		   (*tss).Sq1 = Sq + 1;
		   (*tss).ZGS_col = (*tss).Sq1 * (*tss).stream_snps;
		   (*tss).AF.resize((*tss).stream_snps);


		   (*tss).stream_i = 0;
		   while ((*tss).stream_i < (*tss).stream_snps) {

			       if ((*tss).snploop == (end + 1) && (*tss).stream_i == 0) {
				        break;
			       }

				   if ((*tss).snploop == end + 1 && (*tss).stream_i != 0) {
					   (*tss).stream_snps = (*tss).stream_i;
					   (*tss).Sq1 = Sq + 1;
					   (*tss).ZGS_col = (*tss).Sq1 * (*tss).stream_snps;
					   (*tss).AF.resize((*tss).stream_snps);
					   break;
				   }

			       (*tss).snploop++;


			       // Number of individuals. Only present when Layout == 1
			       if (test.Layout == 1) {
				       fread(&(*tss).Nrow, 4, 1, (*tss).fin3); // cout << "Nrow: " << Nrow << " " << std::flush;  
				       if ((*tss).Nrow != test.Nbgen) {
					        cerr << "ERROR: Nrow = " << (*tss).Nrow << " does not match Nbgen = " << test.Nbgen << '\n';
					        exit(1);
				       }
			       }


				   // The length of the variant identifier
				   fread(&(*tss).LS, 2, 1, (*tss).fin3);  // cout << "LS: " << LS << " " << std::flush;
				   if ((*tss).LS > (*tss).maxLA) {
					   (*tss).maxLA = 2 * (*tss).LS;
					   delete[](*tss).snpID;
					   (*tss).snpID = new char[(*tss).maxLA + 1];
				   }

				   // The variant identifier
				   fread((*tss).snpID, 1, (*tss).LS, (*tss).fin3); (*tss).snpID[(*tss).LS] = '\0'; // cout << "snpID: " << string(snpID) << " " << std::flush;


				   // The length of the rsid
				   fread(&(*tss).LR, 2, 1, (*tss).fin3); // cout << "LR: " << LR << " " << std::flush;
				   if ((*tss).LR > (*tss).maxLA) {
					   (*tss).maxLA = 2 * (*tss).LR;
					   delete[](*tss).rsID;
					   (*tss).rsID = new char[(*tss).maxLA + 1];
				   }

				   // The rsid
				   fread((*tss).rsID, 1, (*tss).LR, (*tss).fin3); (*tss).rsID[(*tss).LR] = '\0'; // cout << "rsID: " << string(rsID) << " " << std::flush;


				   // The length of the chromosome
				   fread(&(*tss).LC, 2, 1, (*tss).fin3); // cout << "LC: " << LC << " " << std::flush;

				   // The chromosome
				   fread((*tss).chrStr, 1, (*tss).LC, (*tss).fin3); (*tss).chrStr[(*tss).LC] = '\0';

				   // The variant position
				   fread(&(*tss).physpos, 4, 1, (*tss).fin3); // cout << "physpos: " << physpos << " " << std::flush;


				   // The number of alleles if Layout == 2. If Layout == 1, this value is assumed to be 2 and not stored
				   if (test.Layout == 2) {
					   fread(&(*tss).LKnum, 2, 1, (*tss).fin3); // this is for Layout = 2, Lnum = 2 is Layout = 1
					   if ((*tss).LKnum != 2) {
						    cerr << "ERROR: Non-bi-allelic variant found: " << (*tss).LKnum << " alleles\n";
						    exit(1);
					   }
				   }

				   // Length of the first allele
				   fread(&(*tss).LA, 4, 1, (*tss).fin3); // cout << "LA: " << LA << " " << std::flush;
				   if ((*tss).LA > (*tss).maxLA) {
                       (*tss).maxLA = 2 * (*tss).LA;
					   delete[](*tss).allele1;
					   (*tss).allele1 = new char[(*tss).maxLA + 1];
				   }
				   // The first allele
				   fread((*tss).allele1, 1, (*tss).LA, (*tss).fin3); (*tss).allele1[(*tss).LA] = '\0';


				   // The length of the second allele
				   fread(&(*tss).LB, 4, 1, (*tss).fin3); // cout << "LB: " << LB << " " << std::flush;
				   if ((*tss).LB > (*tss).maxLB) {
					   (*tss).maxLB = 2 * (*tss).LB;
					   delete[](*tss).allele0;
					   (*tss).allele0 = new char[(*tss).maxLB + 1];
				   }
				   // The second allele
				   fread((*tss).allele0, 1, (*tss).LB, (*tss).fin3); (*tss).allele0[(*tss).LB] = '\0';



				   //cout << string((*tss).snpID) + "\t" + string((*tss).rsID) + "\t" + string((*tss).chrStr) + "\t" + std::to_string((*tss).physpos) + "\t" + string((*tss).allele1) + "\t" + string((*tss).allele0);
				   if (test.Layout == 1) {

					   fread(&(*tss).zLen, 4, 1, (*tss).fin3); // cout << "zLen: " << zLen << endl;
					   fread(&(*tss).zBuf[0], 1, (*tss).zLen, (*tss).fin3);

					   (*tss).destLen = 6 * test.Nbgen;

					   if (uncompress(&(*tss).shortBuf[0], &(*tss).destLen, &(*tss).zBuf[0], (*tss).zLen) != Z_OK || (*tss).destLen != 6 * test.Nbgen) {
						   cerr << "ERROR: uncompress() failed\n";
						   exit(1);
					   }


					   // read genotype probabilities
					   const double scale = 1.0 / 32768;
					   (*tss).tmp1 = (*tss).stream_i * (*tss).Sq1 * samSize;

					   //if (IDMatching == 1) {
					   (*tss).idx_k = 0;
					   for (uint i = 0; i < test.Nbgen; i++) {

						    if (include_idx[(*tss).idx_k] == i) {
							    (*tss).p11 = (*tss).shortBuf[3 * i] * scale;
							    (*tss).p10 = (*tss).shortBuf[3 * i + 1] * scale;
							    (*tss).p00 = (*tss).shortBuf[3 * i + 2] * scale;

							    (*tss).pTot = (*tss).p11 + (*tss).p10 + (*tss).p00;
							    (*tss).dosage = (2 * (*tss).p00 + (*tss).p10) / (*tss).pTot;

							    (*tss).tmp2 = (*tss).idx_k + (*tss).tmp1;
							    (*tss).AF[(*tss).stream_i] += (*tss).dosage;

							    if (test.phenoTyp == 1) {
								   (*tss).ZGSvec[(*tss).tmp2] = test.miu[(*tss).idx_k] * (1 - test.miu[(*tss).idx_k]) * (*tss).dosage;
							    } else {
								   (*tss).ZGSvec[(*tss).tmp2] = (*tss).dosage;
							    }

							   (*tss).idx_k++;
						    }
					   }

					   if (((*tss).AF[(*tss).stream_i] / 2 / samSize) < test.maf || ((*tss).AF[(*tss).stream_i] / 2 / samSize) > (1 - test.maf)) {
						    (*tss).AF[(*tss).stream_i] = 0;
						     continue;
					   }
					   for (int j = 0; j < Sq; j++) {
						    (*tss).tmp3 = samSize * (j + 1) + (*tss).tmp1;

						    for (uint i = 0; i < samSize; i++) {
							     int tmp4 = i * (numSelCol + 1);
								 (*tss).ZGSvec[(*tss).tmp3 + i] = test.covX[tmp4 + j + 1] * (*tss).ZGSvec[(*tss).tmp1 + i]; // here we save ZGS in column wise
						    }

					   }

				   } // end of reading genotype data when Layout = 1




				   if(test.Layout= 2){

					  fread(&(*tss).zLen, 4, 1, (*tss).fin3); // cout << "zLen: " << zLen << endl;
				      (*tss).zBuf.resize((*tss).zLen - 4);

				      if (test.CompressedSNPBlocks == 0) {
					      (*tss).DLen = (*tss).zLen;
					      fread(&(*tss).zBuf[0], 1, (*tss).zLen, (*tss).fin3);
				      } else {
					      fread(&(*tss).DLen, 4, 1, (*tss).fin3);
					      fread(&(*tss).zBuf[0], 1, (*tss).zLen - 4, (*tss).fin3);
				      }


				      (*tss).destLen = (*tss).DLen; //6*Nbgen;
				      (*tss).shortBuf.resize((*tss).DLen);


					  if (uncompress(&(*tss).shortBuf[0], &(*tss).destLen, &(*tss).zBuf[0], (*tss).zLen - 4) != Z_OK || (*tss).destLen != (*tss).DLen) {
						  cout << "destLen: " << (*tss).destLen << " " << (*tss).zLen - 4 << '\n';
						  cerr << "ERROR: uncompress() failed\n";
						  exit(1);
				 	  }


				      // read genotype probabilities
				      (*tss).bufAt = &(*tss).shortBuf[0];
				      (*tss).N = (*tss).bufAt[0] | ((*tss).bufAt[1] << 8) | ((*tss).bufAt[2] << 16) | ((*tss).bufAt[3] << 24); (*tss).bufAt += 4;
				      if ((*tss).N != test.Nbgen) {
					      cerr << "ERROR: " << "snpName " << " has N = " << (*tss).N << " (mismatch with header block)\n";
					      exit(1);
				      }



				      (*tss).K = (*tss).bufAt[0] | ((*tss).bufAt[1] << 8); (*tss).bufAt += 2;
				      if ((*tss).K != 2U) {
					       cout << "\n WARNING: There are SNP(s) with more than 2 alleles (non-bi-allelic). Skipping... \n\n";
exit(1);
					  }


					  (*tss).Pmin = *(*tss).bufAt; (*tss).bufAt++;
					  if ((*tss).Pmin > 2U) {
						  cerr << "ERROR: " << (*tss).snpID << " has minimum ploidy = " << (*tss).Pmin << " (not 2)\n";
						  exit(1);
					  }


					  (*tss).Pmax = *(*tss).bufAt; (*tss).bufAt++;
					  if ((*tss).Pmax != 2U) {
						  cerr << "ERROR: " << (*tss).snpID << " has maximum ploidy = " << (*tss).Pmax << " (not 2)\n";
						  exit(1);
					  }




					  (*tss).ploidy_and_missing_info.resize((*tss).N);
					  (*tss).ploidy_sum = 0;
					  (*tss).idx_to_sum = 0;
					  for (uint i = 0; i < (*tss).N; i++) {
						  (*tss).ploidyMiss = *(*tss).bufAt;
						  (*tss).ploidy_and_missing_info[i] = (*tss).ploidyMiss;

						  if ((*tss).ploidyMiss > 2U) {
							  std::cerr << "ERROR: " << (*tss).snpID << " has ploidy/missingness byte = " << (*tss).ploidyMiss
								  << " (not 2)" << endl;
							  exit(1);
						  }

						  if (test.include_idx[(*tss).idx_to_sum] == i) {
							  (*tss).ploidy_sum += (*tss).ploidyMiss;
							  (*tss).idx_to_sum++;
						  }

						  (*tss).bufAt++;
					  }


					  // Phased information  indicating what is stored in the row. 
					  //    Phased = 1; row stores one probability per allele
					  //    Phased = 0; row stores one probability per possible genotype
					  //    Everything else is error.
					  (*tss).Phased = *(*tss).bufAt; (*tss).bufAt++;
					  if ((*tss).Phased != 0U && (*tss).Phased != 1U) {
						  //if(Phased == 1U){ cerr << "\nERROR: Phased data is not currently supported by GEM.\n"; exit(1);}
						  cerr << "\nERROR: " << (*tss).snpID << " has Phased = " << (*tss).Phased << ". Must be 0 or 1. \n"
							  << "       See https://www.well.ox.ac.uk/~gav/bgen_format/spec/latest.html for more details." << endl;
						  exit(1);
					  }

					  (*tss).B = *(*tss).bufAt; (*tss).bufAt++;
					  (*tss).Bbits = std::pow(2, (*tss).B);
					  if (((*tss).B != 8U) && ((*tss).B != 16U) && ((*tss).B != 24U) && ((*tss).B != 32U)) {
						  std::cerr << "ERROR: " << "snpName " << " has B = " << (*tss).B << " (not divisible by 8)\n";
						  exit(1);
					  }



					  (*tss).tmp1 = (*tss).stream_i * (*tss).Sq1 * samSize;
					  (*tss).idx_k = 0;

					  for (uint i = 0; i < (*tss).N; i++) {

						  if ((*tss).B == 8U)
							  (*tss).chartem = (*tss).bufAt[0];
						  else if ((*tss).B == 16U)
							  (*tss).chartem = (*tss).bufAt[0] | ((*tss).bufAt[1] << 8);
						  else if ((*tss).B == 24U)
							  (*tss).chartem = (*tss).bufAt[0] | ((*tss).bufAt[1] << 8) | ((*tss).bufAt[2] << 16);
						  else if ((*tss).B == 32U)
							  (*tss).chartem = (*tss).bufAt[0] | ((*tss).bufAt[1] << 8) | ((*tss).bufAt[2] << 16) | ((*tss).bufAt[3] << 24);
						  (*tss).bufAt += (*tss).B / 8;


						  if ((*tss).ploidy_and_missing_info[i] == 2U) {
							  if ((*tss).B == 8U)
								  (*tss).chartem1 = (*tss).bufAt[0];
							  else if ((*tss).B == 16U)
								  (*tss).chartem1 = (*tss).bufAt[0] | ((*tss).bufAt[1] << 8);
							  else if ((*tss).B == 24U)
								  (*tss).chartem1 = (*tss).bufAt[0] | ((*tss).bufAt[1] << 8) | ((*tss).bufAt[2] << 16);
							  else if ((*tss).B == 32U)
								  (*tss).chartem1 = (*tss).bufAt[0] | ((*tss).bufAt[1] << 8) | ((*tss).bufAt[2] << 16) | ((*tss).bufAt[3] << 24);
							  (*tss).bufAt += (*tss).B / 8;
						  }




						  if (include_idx[(*tss).idx_k] == i) {
							  (*tss).p11 = (*tss).chartem / double(1.0 * ((*tss).Bbits - 1));
							  (*tss).p10 = (*tss).chartem1 / double(1.0 * ((*tss).Bbits - 1));
							  (*tss).dosage;

							  if ((*tss).Phased == 1U) {
								   if((*tss).ploidy_and_missing_info[i] == 1U){
									   (*tss).dosage = 1 - (*tss).p11;
								   } else {
									   (*tss).dosage = 2 - ((*tss).p11 + (*tss).p10);
								   }

							   } else {
								   (*tss).dosage = 2 * (1 - (*tss).p11 - (*tss).p10) + (*tss).p10;
							   }

							   (*tss).tmp2 = (*tss).idx_k + (*tss).tmp1;
							   (*tss).AF[(*tss).stream_i] += (*tss).dosage;

							   if (test.phenoTyp == 1) {
								   (*tss).ZGSvec[(*tss).tmp2] = test.miu[(*tss).idx_k] * (1 - test.miu[(*tss).idx_k]) * (*tss).dosage;
							   } else {
                                   (*tss).ZGSvec[(*tss).tmp2] = (*tss).dosage; // replace your new data from other genotype files here.
							   }

							   (*tss).idx_k++;
						   }

				      }

				      // cout << string((*tss).snpID) + "\t" + string((*tss).rsID) + "\t" + string((*tss).chrStr) + "\t" + std::to_string((*tss).physpos) + "\t" + string((*tss).allele1) + "\t" + string((*tss).allele0) << "\t AF: " << (*tss).AF[(*tss).stream_i] << endl;
				      if (((*tss).AF[(*tss).stream_i] / 2 / samSize) < test.maf || ((*tss).AF[(*tss).stream_i] / 2 / samSize) > (1 - test.maf)) {
					      (*tss).AF[(*tss).stream_i] = 0;
					      continue;
				      }

				      (*tss).geno_snpid[(*tss).stream_i] = string((*tss).snpID) + "\t" + string((*tss).rsID) + "\t" + string((*tss).chrStr) + "\t" + std::to_string((*tss).physpos) + "\t" + string((*tss).allele1) + "\t" + string((*tss).allele0) + "\t" + std::to_string((*tss).AF[(*tss).stream_i] / (*tss).ploidy_sum);


					  for (int j = 0; j < Sq; j++) {
						   (*tss).tmp3 = samSize * (j + 1) + (*tss).tmp1;
						   for (uint i = 0; i < samSize; i++) {
							    (*tss).tmp4 = i * (numSelCol + 1);
							    (*tss).ZGSvec[(*tss).tmp3 + i] = test.covX[(*tss).tmp4 + j + 1] * (*tss).ZGSvec[(*tss).tmp1 + i]; // here we save ZGS in column wise
						   }
					  }


			       } // end of layout 2

				   (*tss).stream_i++;
		   } // end of stream_i

		   if ((*tss).snploop == end + 1 & (*tss).stream_i == 0) {
			   break;
		   }


		   /***************************************************************/
		   //	genodata and envirment data
		   (*tss).ZGS = &(*tss).ZGSvec[0];

		   // transpose(X) * ZGS
		   // it is non-squred matrix, attention that continuous memory is column-major due to Fortran in BLAS.
		   // important!!!!
		   (*tss).XtransZGS = new double[(numSelCol + 1) * (*tss).ZGS_col];
		   matNmatNprod(test.covX, (*tss).ZGS, (*tss).XtransZGS, numSelCol + 1, samSize, (*tss).ZGS_col);



		   if (test.phenoTyp == 0) {
			   for (int j = 0; j < (*tss).ZGS_col; j++) {
				    for (int i = 0; i < samSize; i++) {
					     for (int k = 0; k <= numSelCol; k++) {
						     (*tss).ZGS[j * samSize + i] -= test.XinvXTX[k * samSize + i] * (*tss).XtransZGS[j * (numSelCol + 1) + k];
					     }
					     if (robust == 1) { (*tss).ZGSR2vec[j * samSize + i] = (*tss).ZGS[j * samSize + i] * test.resid[i] * test.resid[i]; }
				    }
			   }
		   }
		   else if (test.phenoTyp == 1) {

			   for (int j = 0; j < (*tss).ZGS_col; j++) {
					for (int i = 0; i < samSize; i++) {
						(*tss).ZGStemp = 0.0;
						for (int k = 0; k <= numSelCol; k++) {
							if (k == 0) {
								(*tss).ZGStemp += (*tss).ZGS[j * samSize + i] / test.miu[i] / (1.0 - test.miu[i]) - test.XinvXTX[k * samSize + i] * (*tss).XtransZGS[j * (numSelCol + 1) + k];
							}
							else {
								(*tss).ZGStemp -= test.XinvXTX[k * samSize + i] * (*tss).XtransZGS[j * (numSelCol + 1) + k];
							}
						}

						(*tss).ZGS[j * samSize + i] = (*tss).ZGStemp;
						if (robust == 1) (*tss).ZGSR2vec[j * samSize + i] = (*tss).ZGS[j * samSize + i] * test.resid[i] * test.resid[i];
						(*tss).WZGS[j * samSize + i] = test.miu[i] * (1 - test.miu[i]) * (*tss).ZGS[j * samSize + i];
					}
			   }
		   }





		   (*tss).ZGSR2 = &(*tss).ZGSR2vec[0];
		   delete[](*tss).XtransZGS;
		   // transpose(ZGS) * resid
		   (*tss).ZGStR = new double[(*tss).ZGS_col];
		   matvecprod((*tss).ZGS, test.resid, (*tss).ZGStR, (*tss).ZGS_col, samSize);
		   // transpose(ZGS) * ZGS
		   (*tss).ZGStZGS = new double[(*tss).ZGS_col * (*tss).ZGS_col];
		   if (test.phenoTyp == 0) {
			   matmatTprod((*tss).ZGS, (*tss).ZGS, (*tss).ZGStZGS, (*tss).ZGS_col, samSize, (*tss).ZGS_col);
		   }
		   else if (test.phenoTyp == 1) {
			   matmatTprod((*tss).ZGS, (*tss).WZGS, (*tss).ZGStZGS, (*tss).ZGS_col, samSize, (*tss).ZGS_col);
		   }
		   else {
			   cout << "phenoTyp is not equal to 0 or 1. Kill the job!! \n";
			   exit(1);
		   }


		   // transpose(ZGSR2) * ZGS
		   (*tss).ZGSR2tZGS = new double[(*tss).ZGS_col * (*tss).ZGS_col];
		   if (robust == 1) matmatTprod((*tss).ZGSR2, (*tss).ZGS, (*tss).ZGSR2tZGS, (*tss).ZGS_col, samSize, (*tss).ZGS_col);

		   (*tss).betaM = new double[(*tss).stream_snps];
		   (*tss).VarbetaM = new double[(*tss).stream_snps];
		   (*tss).betaInt = new double* [(*tss).stream_snps];
		   (*tss).VarbetaInt = new double* [(*tss).stream_snps];
		   (*tss).PvalM = new double[(*tss).stream_snps];
		   (*tss).PvalInt = new double[(*tss).stream_snps];
		   (*tss).PvalJoint = new double[(*tss).stream_snps];
		   boost::math::chi_squared chisq_dist_M(1);
		   boost::math::chi_squared chisq_dist_Int(Sq);
		   boost::math::chi_squared chisq_dist_Joint(1 + Sq);


		   if (robust == 0) {

				for (int i = 0; i < (*tss).stream_snps; i++) {
					// initialize dynamic 2D array
					(*tss).betaInt[i] = new double[Sq];
					(*tss).VarbetaInt[i] = new double[Sq * Sq];

					// betamain
					(*tss).tmp1 = i * (*tss).ZGS_col * (*tss).Sq1 + i * (*tss).Sq1;
					(*tss).betaM[i] = (*tss).ZGStR[i * (*tss).Sq1] / (*tss).ZGStZGS[(*tss).tmp1];
					(*tss).VarbetaM[i] = test.sigma2 / (*tss).ZGStZGS[(*tss).tmp1];

					(*tss).S2TransS2 = new double[Sq * Sq];
					(*tss).S2TransR = new double[Sq];
					(*tss).S2DS2 = new double[Sq * Sq];
					(*tss).InvVarbetaint = new double[Sq * Sq];
					for (int ind1 = 0; ind1 < Sq; ind1++) {
						for (int ind2 = 0; ind2 < Sq; ind2++) {
							// transpose(Snew2) * Snew2
							(*tss).S2TransS2[ind1 * Sq + ind2] = (*tss).ZGStZGS[(*tss).tmp1 + (ind1 + 1) * (*tss).ZGS_col + ind2 + 1] - (*tss).ZGStZGS[(*tss).tmp1 + (ind1 + 1) * (*tss).ZGS_col] * (*tss).ZGStZGS[(*tss).tmp1 + ind2 + 1] / (*tss).ZGStZGS[(*tss).tmp1];
						}
						//transpose(Snew2) * resid
						(*tss).S2TransR[ind1] = (*tss).ZGStR[i * (*tss).Sq1 + ind1 + 1] - (*tss).ZGStZGS[(*tss).tmp1 + (ind1 + 1) * (*tss).ZGS_col] * (*tss).ZGStR[i * (*tss).Sq1] / (*tss).ZGStZGS[(*tss).tmp1];
					}

					// invert (S2TransS2)
					matInv((*tss).S2TransS2, Sq);

					// betaInt = invert(S2TransS2) * S2TransR
					matvecprod((*tss).S2TransS2, (*tss).S2TransR, (*tss).betaInt[i], Sq, Sq);

					// Inv(S2TransS2) * S2DS2
					(*tss).Stemp2 = new double[Sq * Sq];

					for (int j = 0; j < Sq; j++) {
						for (int k = 0; k < Sq; k++) {
							(*tss).VarbetaInt[i][j * Sq + k] = test.sigma2 * (*tss).S2TransS2[j * Sq + k];
							(*tss).InvVarbetaint[j * Sq + k] = (*tss).VarbetaInt[i][j * Sq + k];
						}
					}

					// calculating P values
					(*tss).statM = (*tss).betaM[i] * (*tss).betaM[i] / (*tss).VarbetaM[i];
					if (isnan((*tss).statM) || (*tss).statM <= 0.0) {
						(*tss).PvalM[i] = NAN;
					}
					else {
						(*tss).PvalM[i] = boost::math::cdf(complement(chisq_dist_M, (*tss).statM));
					}

					// invert VarbetaInt[i]
					matInv((*tss).InvVarbetaint, Sq);
					(*tss).Stemp3 = new double[Sq];
					matvecprod((*tss).InvVarbetaint, (*tss).betaInt[i], (*tss).Stemp3, Sq, Sq);

					(*tss).statInt = 0.0;
					for (int j = 0; j < Sq; j++) {
						(*tss).statInt += (*tss).betaInt[i][j] * (*tss).Stemp3[j];
					}

					if (isnan((*tss).statInt) || (*tss).statInt <= 0.0) {
						(*tss).PvalInt[i] = NAN;
					}
					else {
						(*tss).PvalInt[i] = boost::math::cdf(complement(chisq_dist_Int, (*tss).statInt));
					}


					(*tss).statJoint = (*tss).statM + (*tss).statInt;
					if (isnan((*tss).statJoint) || (*tss).statJoint <= 0.0) {
						(*tss).PvalJoint[i] = NAN;
					}
					else {
						(*tss).PvalJoint[i] = boost::math::cdf(complement(chisq_dist_Joint, (*tss).statJoint));
					}

					delete[](*tss).S2TransS2;
					delete[](*tss).S2TransR;
					delete[](*tss).S2DS2;
					delete[](*tss).Stemp2;
					delete[](*tss).Stemp3;
					delete[](*tss).InvVarbetaint;
				}
		   }
		   else if (robust == 1) {

				for (int i = 0; i < (*tss).stream_snps; i++) {
					// initialize dynamic 2D array
					(*tss).betaInt[i] = new double[Sq];
					(*tss).VarbetaInt[i] = new double[Sq * Sq];

					//betamain
					(*tss).tmp1 = i * (*tss).ZGS_col * (*tss).Sq1 + i * (*tss).Sq1;
					(*tss).betaM[i] = (*tss).ZGStR[i * (*tss).Sq1] / (*tss).ZGStZGS[(*tss).tmp1];
					(*tss).VarbetaM[i] = (*tss).ZGSR2tZGS[(*tss).tmp1] / ((*tss).ZGStZGS[(*tss).tmp1] * (*tss).ZGStZGS[(*tss).tmp1]);

					(*tss).S2TransS2 = new double[Sq * Sq];
					(*tss).S2TransR = new double[Sq];
					(*tss).S2DS2 = new double[Sq * Sq];
					(*tss).InvVarbetaint = new double[Sq * Sq];

					for (int ind1 = 0; ind1 < Sq; ind1++) {
						for (int ind2 = 0; ind2 < Sq; ind2++) {
							// transpose(Snew2) * Snew2
							(*tss).S2TransS2[ind1 * Sq + ind2] = (*tss).ZGStZGS[(*tss).tmp1 + (ind1 + 1) * (*tss).ZGS_col + ind2 + 1] - (*tss).ZGStZGS[(*tss).tmp1 + (ind1 + 1) * (*tss).ZGS_col] * (*tss).ZGStZGS[(*tss).tmp1 + ind2 + 1] / (*tss).ZGStZGS[(*tss).tmp1];
							// transpose(Snew2) * D * Snew2
							(*tss).S2DS2[ind1 * Sq + ind2] = (*tss).ZGSR2tZGS[(*tss).tmp1 + (ind1 + 1) * (*tss).ZGS_col + ind2 + 1] - (*tss).ZGStZGS[(*tss).tmp1 + (ind1 + 1) * (*tss).ZGS_col] * (*tss).ZGSR2tZGS[(*tss).tmp1 + ind2 + 1] / (*tss).ZGStZGS[(*tss).tmp1] - (*tss).ZGSR2tZGS[(*tss).tmp1 + (ind1 + 1) * (*tss).ZGS_col] * (*tss).ZGStZGS[(*tss).tmp1 + ind2 + 1] / (*tss).ZGStZGS[(*tss).tmp1] + (*tss).ZGStZGS[(*tss).tmp1 + (ind1 + 1) * (*tss).ZGS_col] * (*tss).ZGSR2tZGS[(*tss).tmp1] * (*tss).ZGStZGS[(*tss).tmp1 + ind2 + 1] / (*tss).ZGStZGS[(*tss).tmp1] / (*tss).ZGStZGS[(*tss).tmp1];
						}
						//transpose(Snew2) * resid
						(*tss).S2TransR[ind1] = (*tss).ZGStR[i * (*tss).Sq1 + ind1 + 1] - (*tss).ZGStZGS[(*tss).tmp1 + (ind1 + 1) * (*tss).ZGS_col] * (*tss).ZGStR[i * (*tss).Sq1] / (*tss).ZGStZGS[(*tss).tmp1];
					}

					// invert (S2TransS2)
					matInv((*tss).S2TransS2, Sq);

					// betaInt = invert(S2TransS2) * S2TransR
					matvecprod((*tss).S2TransS2, (*tss).S2TransR, (*tss).betaInt[i], Sq, Sq);

					// Inv(S2TransS2) * S2DS2
					(*tss).Stemp2 = new double[Sq * Sq];
					matmatprod((*tss).S2TransS2, (*tss).S2DS2, (*tss).Stemp2, Sq, Sq, Sq);

					// Stemp2 * Inv(S2TransS2)
					matNmatNprod((*tss).Stemp2, (*tss).S2TransS2, (*tss).VarbetaInt[i], Sq, Sq, Sq);



					for (int j = 0; j < Sq; j++) {
						for (int k = 0; k < Sq; k++) {
							(*tss).InvVarbetaint[j * Sq + k] = (*tss).VarbetaInt[i][j * Sq + k];
						}
					}


					//calculating P values
					(*tss).statM = (*tss).betaM[i] * (*tss).betaM[i] / (*tss).VarbetaM[i];
					if (isnan((*tss).statM) || (*tss).statM <= 0.0) {
						(*tss).PvalM[i] = NAN;
					}
					else {
						(*tss).PvalM[i] = boost::math::cdf(complement(chisq_dist_M, (*tss).statM));
					}


					// invert VarbetaInt[i]
					matInv((*tss).InvVarbetaint, Sq);
					(*tss).Stemp3 = new double[Sq];
					matvecprod((*tss).InvVarbetaint, (*tss).betaInt[i], (*tss).Stemp3, Sq, Sq);

					(*tss).statInt = 0.0;
					for (int j = 0; j < Sq; j++) {
						(*tss).statInt += (*tss).betaInt[i][j] * (*tss).Stemp3[j];
					}

					if (isnan((*tss).statInt) || (*tss).statInt <= 0.0) {
						(*tss).PvalInt[i] = NAN;
					}
					else {
						(*tss).PvalInt[i] = boost::math::cdf(complement(chisq_dist_Int, (*tss).statInt));
					}

					(*tss).statJoint = (*tss).statM + (*tss).statInt;
					if (isnan((*tss).statJoint) || (*tss).statJoint <= 0.0) {
						(*tss).PvalJoint[i] = NAN;
					}
					else {
						(*tss).PvalJoint[i] = boost::math::cdf(complement(chisq_dist_Joint, (*tss).statJoint));
					}


					delete[](*tss).S2TransS2;
					delete[](*tss).S2TransR;
					delete[](*tss).S2DS2;
					delete[](*tss).Stemp2;
					delete[](*tss).Stemp3;
					delete[](*tss).InvVarbetaint;
				}
		   } // end of if robust == 1

		  

		   for (int i = 0; i < (*tss).stream_snps; i++) {
			    results << (*tss).geno_snpid[i] << "\t" << (*tss).betaM[i] << "\t" << (*tss).VarbetaM[i] << "\t";
				for (int ii = 0; ii < Sq; ii++) {
					results << (*tss).betaInt[i][ii] << "\t";
				}

				for (int ii = 0; ii < Sq; ii++) {
					for (int jj = 0; jj < Sq; jj++) {
						results << (*tss).VarbetaInt[i][ii * Sq + jj] << "\t";
					}
				}
			    results << (*tss).PvalM[i] << "\t" << (*tss).PvalInt[i] << "\t" << (*tss).PvalJoint[i] << '\n';
		   }

		   delete[](*tss).ZGStR;
		   delete[](*tss).ZGStZGS;
		   delete[](*tss).ZGSR2tZGS;

		   delete[](*tss).betaM;
		   delete[](*tss).VarbetaM;
		   for (int i = 0; i < (*tss).stream_snps; i++) {
			    delete[](*tss).betaInt[i];
			    delete[](*tss).VarbetaInt[i];
		   }
		   delete[](*tss).betaInt;
		   delete[](*tss).VarbetaInt;
		   delete[](*tss).PvalM;
		   delete[](*tss).PvalInt;
		   delete[](*tss).PvalJoint;
		   (*tss).AF.clear();
	} // end of snploop 


	// Close files
	results.close();
	fclose((*tss).fin3);

}