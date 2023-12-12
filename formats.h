#ifndef FORMATS_H_
#define FORMATS_H_
#include <iostream>
enum file_format
{
	FASTA = 1,
	FASTA_CONT,
	FASTQ,
	BAM,
	TAB_MATE5,
	TAB_MATE6,
	RAW,
	CMDLINE,
	QSEQ,
	SRA_FASTA,
	SRA_FASTQ
};
static const std::string file_format_names[] = {
	"Invalid!",
	"FASTA",
	"FASTA sampling",
	"FASTQ",
	"BAM",
	"Tabbed 5-field",
	"Tabbed 6-field",
	"Raw",
	"Command line",
	"Qseq",
	"SRA Fasta",
	"SRA Fastq"};
#endif
