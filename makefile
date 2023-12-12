gFLAG := -O3 -g -Wliteral-suffix -Wall -DNDEBUG -DPOPCNT_CAPABILITY -std=c++11 -msse2
effaln: aln_search.o main.o qual.o pat.o sam.o \
  read_qseq.o aligner_seed_policy.o aligner_seed.o aligner_seed2.o aligner_sw.o \
  aligner_sw_driver.o aligner_cache.o aligner_result.o ref_coord.o mask.o \
  pe.o aln_sink.o dp_framer.o scoring.o presets.o unique.o simple_func.o random_util.o \
  aligner_bt.o sse_util.o aligner_swsse.o outq.o \
  aligner_driver.o tinythread.o ccnt_lut.o ref_read.o alphabet.o shmem.o edit.o aln_idx.o aln_io.o aln_util.o \
  reference.o ds.o multikey_qsort.o limit.o random_source.o embedding.o InArray.o BitMap.o loadkit.o savekit.o
	g++ $(gFLAG) -L ./semiWFA/build -o $@ $^ -lpthread -lz -no-pie -lwfa -ltbb
effaln-index: aln_build.o build_main.o qual.o pat.o sam.o read_qseq.o \
  aligner_cache.o aligner_result.o ref_coord.o mask.o diff_sample.o\
  pe.o scoring.o presets.o unique.o simple_func.o random_util.o sse_util.o \
  tinythread.o ccnt_lut.o ref_read.o alphabet.o shmem.o edit.o aln_idx.o aln_io.o aln_util.o \
  reference.o ds.o multikey_qsort.o limit.o random_source.o InArray.o BitMap.o loadkit.o savekit.o
	g++ $(gFLAG) -o $@ $^ -lpthread -lz -no-pie -ltbb
aln_build.o: aln_build.cpp aln_idx.h alphabet.h assert_helpers.h endian_swap.h aln_idx.h formats.h sequence_io.h tokenize.h aln_sink.h pat.h threading.h ds.h aligner_metrics.h sam.h aligner_seed.h aligner_seed_policy.h aligner_driver.h aligner_sw.h aligner_sw_driver.h aligner_cache.h util.h opts.h outq.h aligner_seed2.h
	g++ $(gFLAG) -I./ -DWITH_TBB -DBUILD_HOST="\"$${HOSTNAME:-`hostname`}\"" -DBUILD_TIME="\"`date -u`\"" -DCOMPILER_VERSION="\"`g++ -v 2>&1 | tail -1`\"" -DCOMPILER_OPTIONS -lpthread -ltbb -lz -c $<
diff_sample.o: diff_sample.cpp diff_sample.h
	g++ $(gFLAG) -I./ -c $<
build_main.o: build_main.cpp tokenize.h ds.h mem_ids.h
	g++ $(gFLAG) -c $<
aln_search.o: aln_search.cpp aln_search.h alphabet.h assert_helpers.h endian_swap.h aln_idx.h formats.h sequence_io.h tokenize.h aln_sink.h pat.h threading.h ds.h aligner_metrics.h sam.h aligner_seed.h aligner_seed_policy.h aligner_driver.h aligner_sw.h aligner_sw_driver.h aligner_cache.h util.h opts.h outq.h aligner_seed2.h aln_search.h
	g++ $(gFLAG) -I./ -DWITH_TBB -DBUILD_HOST="\"$${HOSTNAME:-`hostname`}\"" -DBUILD_TIME="\"`date -u`\"" -DCOMPILER_VERSION="\"`g++ -v 2>&1 | tail -1`\"" -DCOMPILER_OPTIONS -lpthread -ltbb -lz -c $<
aligner_sw.o: aligner_sw.cpp aligner_sw.h aligner_result.h search_globals.h scoring.h mask.h semiWFA/gap_affine/affine_wavefront_align.h embedding.h
	g++ $(gFLAG) -I./ -c $<
main.o: main.cpp tokenize.h ds.h
	g++ $(gFLAG) -c $<
aligner_driver.o: aligner_driver.cpp aligner_driver.h
	g++ $(gFLAG) -c $<
aligner_sw_driver.o: aligner_sw_driver.cpp aligner_sw_driver.h
	g++ $(gFLAG) -I./ -c $<
aln_idx.o: aln_idx.cpp aln_idx.h
	g++ $(gFLAG) -c $<
aln_io.o: aln_io.cpp aln_idx.h
	g++ $(gFLAG) -c $<
aln_util.o: aln_util.cpp aln_idx.h
	g++ $(gFLAG) -c $<
embedding.o: embedding.cpp embedding.h
	g++ $(gFLAG) -I./ -c $<
aligner_cache.o: aligner_cache.cpp aligner_cache.h
	g++ $(gFLAG) -c $<
aligner_result.o: aligner_result.cpp aligner_result.h
	g++ $(gFLAG) -c $<
ref_coord.o: ref_coord.cpp ref_coord.h
	g++ $(gFLAG) -c $<
aligner_seed_policy.o: aligner_seed_policy.cpp aligner_seed_policy.h
	g++ $(gFLAG) -c $<
aligner_seed.o: aligner_seed.cpp aligner_seed.h
	g++ $(gFLAG) -c $<
aligner_seed2.o: aligner_seed2.cpp aligner_seed2.h
	g++ $(gFLAG) -I./ -c $<
mask.o: mask.cpp mask.h
	g++ $(gFLAG) -c $<
pe.o: pe.cpp pe.h
	g++ $(gFLAG) -c $<
aln_sink.o: aln_sink.cpp aln_sink.h
	g++ $(gFLAG) -c $<
dp_framer.o: dp_framer.cpp dp_framer.h
	g++ $(gFLAG) -c $<
scoring.o: scoring.cpp scoring.h
	g++ $(gFLAG) -c $<
presets.o: presets.cpp presets.h opts.h ds.h
	g++ $(gFLAG) -c $<
unique.o: unique.cpp unique.h
	g++ $(gFLAG) -c $<
simple_func.o: simple_func.cpp simple_func.h
	g++ $(gFLAG) -c $<
random_util.o: random_util.cpp random_util.h
	g++ $(gFLAG) -c $<
aligner_bt.o: aligner_bt.cpp aligner_bt.h
	g++ $(gFLAG) -c $<
sse_util.o: sse_util.cpp sse_util.h
	g++ $(gFLAG) -c $<
aligner_swsse.o: aligner_swsse.cpp aligner_swsse.h
	g++ $(gFLAG) -c $<
outq.o: outq.cpp outq.h
	g++ $(gFLAG) -c $<
qual.o: qual.cpp qual.h
	g++ $(gFLAG) -c $<
pat.o: pat.cpp pat.h
	g++ $(gFLAG) -c $<
sam.o: sam.cpp sam.h
	g++ $(gFLAG) -c $<
read_qseq.o: read_qseq.cpp pat.h
	g++ $(gFLAG) -c $<
tinythread.o: tinythread.cpp tinythread.h fast_mutex.h
	g++ $(gFLAG) -c $<
ccnt_lut.o: ccnt_lut.cpp
	g++ $(gFLAG) -c $<
ref_read.o: ref_read.cpp ref_read.h
	g++ $(gFLAG) -c $<
alphabet.o: alphabet.cpp alphabet.h
	g++ $(gFLAG) -c $<
shmem.o: shmem.cpp shmem.h
	g++ $(gFLAG) -c $<
edit.o: edit.cpp edit.h
	g++ $(gFLAG) -c $<
reference.o: reference.cpp reference.h mem_ids.h
	g++ $(gFLAG) -c $<
ds.o: ds.cpp ds.h ref_coord.h
	g++ $(gFLAG) -c $<
multikey_qsort.o: multikey_qsort.cpp multikey_qsort.h
	g++ $(gFLAG) -c $<
limit.o: limit.cpp limit.h
	g++ $(gFLAG) -c $<
random_source.o: random_source.cpp random_source.h
	g++ $(gFLAG) -c $<
InArray.o: InArray.cpp InArray.h
	g++ $(gFLAG) -c $<
BitMap.o: BitMap.cpp BitMap.h
	g++ $(gFLAG) -c $<
loadkit.o: loadkit.cpp loadkit.h
	g++ $(gFLAG) -c $<
savekit.o: savekit.cpp savekit.h
	g++ $(gFLAG) -c $<
clean:
	rm -rf effaln
cleanall:
	rm -rf effaln *.o