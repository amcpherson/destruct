###################################################
#
# Makefile for tools
#
###################################################

#
# Macros
#

ZDIR = ../external/zlib-1.2.5
ZLIB = libz.a

BAMTOOLSDIR = ../external/BamTools-2.0/

INCLUDE = -I. -I../external -I$(ZDIR) -I$(BAMTOOLSDIR)
CC = g++
CXX = g++
CXXFLAGS = -O3 -ggdb -DDEBUG_CHECKS $(INCLUDE)
LDFLAGS = -O3 -ggdb -L$(ZDIR) -lz -lboost_serialization -lboost_iostreams

BINS = \
		overlapclusters\
		mclustermatepairs\
		mpredictbreaks\
		clustermatepairs\
		predictbreaks\
		presetcover\
		setcover\
		setcover_dnarna\
		initsplitalign\
		dosplitalign\
		evalsplitalign\
		splitseq\
		dosplitaligngen\
		decisionmetrics\
		calcspandist\
		calccov\
		estislands\
		localalign\
		realign2\
		realign\
		realign_rna\
		realignp\
		aligntrue\
		alignnull\
		rankclusters\
		cycles\
		paths\
		filterreads\
		filterrnaseq\
		clusterconnectivity\
		qualtrimfastq\
		testssealign\
		samplefastq\
		bamdiscordantfastq\
		bammixture\
		bamsamplespanning\
		bamsamplefastq\
		bamfastq\
		bamconcordantcounts\
		bamconcordantreads\
		bamstats\
		samplegc\
		samplemap\
		estimategc\
		samplemappable

BAMTOOLS_OBJS = \
		$(BAMTOOLSDIR)/api/BamAlignment.o\
		$(BAMTOOLSDIR)/api/BamMultiReader.o\
		$(BAMTOOLSDIR)/api/BamReader.o\
		$(BAMTOOLSDIR)/api/BamWriter.o\
		$(BAMTOOLSDIR)/api/SamHeader.o\
		$(BAMTOOLSDIR)/api/SamProgram.o\
		$(BAMTOOLSDIR)/api/SamProgramChain.o\
		$(BAMTOOLSDIR)/api/SamReadGroup.o\
		$(BAMTOOLSDIR)/api/SamReadGroupDictionary.o\
		$(BAMTOOLSDIR)/api/SamSequence.o\
		$(BAMTOOLSDIR)/api/SamSequenceDictionary.o\
		$(BAMTOOLSDIR)/api/internal/bam/BamHeader_p.o\
		$(BAMTOOLSDIR)/api/internal/bam/BamMultiReader_p.o\
		$(BAMTOOLSDIR)/api/internal/bam/BamRandomAccessController_p.o\
		$(BAMTOOLSDIR)/api/internal/bam/BamReader_p.o\
		$(BAMTOOLSDIR)/api/internal/bam/BamWriter_p.o\
		$(BAMTOOLSDIR)/api/internal/index/BamIndexFactory_p.o\
		$(BAMTOOLSDIR)/api/internal/index/BamStandardIndex_p.o\
		$(BAMTOOLSDIR)/api/internal/index/BamToolsIndex_p.o\
		$(BAMTOOLSDIR)/api/internal/io/BamDeviceFactory_p.o\
		$(BAMTOOLSDIR)/api/internal/io/BamFile_p.o\
		$(BAMTOOLSDIR)/api/internal/io/BamFtp_p.o\
		$(BAMTOOLSDIR)/api/internal/io/BamHttp_p.o\
		$(BAMTOOLSDIR)/api/internal/io/BamPipe_p.o\
		$(BAMTOOLSDIR)/api/internal/io/BgzfStream_p.o\
		$(BAMTOOLSDIR)/api/internal/io/ByteArray_p.o\
		$(BAMTOOLSDIR)/api/internal/io/HostAddress_p.o\
		$(BAMTOOLSDIR)/api/internal/io/HostInfo_p.o\
		$(BAMTOOLSDIR)/api/internal/io/HttpHeader_p.o\
		$(BAMTOOLSDIR)/api/internal/io/ILocalIODevice_p.o\
		$(BAMTOOLSDIR)/api/internal/io/RollingBuffer_p.o\
		$(BAMTOOLSDIR)/api/internal/io/TcpSocket_p.o\
		$(BAMTOOLSDIR)/api/internal/io/TcpSocketEngine_p.o\
		$(BAMTOOLSDIR)/api/internal/io/TcpSocketEngine_unix_p.o\
		$(BAMTOOLSDIR)/api/internal/sam/SamFormatParser_p.o\
		$(BAMTOOLSDIR)/api/internal/sam/SamFormatPrinter_p.o\
		$(BAMTOOLSDIR)/api/internal/sam/SamHeaderValidator_p.o\
		$(BAMTOOLSDIR)/api/internal/utils/BamException_p.o\
		$(BAMTOOLSDIR)/utils/bamtools_pileup_engine.o

all : $(ZDIR)/$(ZLIB) $(BINS)

#
# Build libs
#

$(ZDIR)/$(ZLIB) :
		$(MAKE) --no-print-directory -C $(ZDIR) $(ZLIB)

#
# Build tools
#

overlapclusters : \
		./Common.o\
		./CompactBreakRegion.o\
		./DebugCheck.o\
		./overlapclusters.o

mclustermatepairs : \
		./AlignmentStream.o\
		./Common.o\
		./DebugCheck.o\
		./MatePairGibbs.o\
		./Parsers.o\
		./asa136.o\
		./asa241.o\
		../external/Triangle/triangle.o\
		./mclustermatepairs.o

mpredictbreaks : \
		./AlignmentStream.o\
		./Common.o\
		./DebugCheck.o\
		./Parsers.o\
		./Sequences.o\
		./mpredictbreaks.o

clustermatepairs : \
		./AlignmentStream.o\
		./Common.o\
		./DebugCheck.o\
		./MatePairGibbs.o\
		./Parsers.o\
		./asa136.o\
		./asa241.o\
		./clustermatepairs.o\
		../external/Triangle/triangle.o

predictbreaks : \
		./AlignmentStream.o\
		./Common.o\
		./DebugCheck.o\
		./Parsers.o\
		./Sequences.o\
		./predictbreaks.o

presetcover : \
		./Common.o\
		./DebugCheck.o\
		./Parsers.o\
		./presetcover.o

setcover : \
  		./Algorithms.o\
		./Common.o\
		./DebugCheck.o\
		./Parsers.o\
		./setcover.o

setcover_dnarna : \
  		./Algorithms.o\
		./Common.o\
		./DebugCheck.o\
		./Parsers.o\
		./setcover_dnarna.o

initsplitalign : \
		./AlignmentStream.o\
		./Common.o\
		./DebugCheck.o\
		./ExonRegions.o\
		./ReadIndex.o\
		./Sequences.o\
		./SplitAlignment.o\
		./SplitReadAligner.o\
		./Parsers.o\
		./initsplitalign.o

dosplitalign : \
		./AlignmentStream.o\
		./ReadIndex.o\
		./Common.o\
		./DebugCheck.o\
		./ExonRegions.o\
		./ReadIndex.o\
		./Sequences.o\
		./SplitAlignment.o\
		./SplitReadAligner.o\
		./Parsers.o\
		./ReadStream.o\
		./dosplitalign.o

evalsplitalign : \
		./AlignmentStream.o\
		./ReadIndex.o\
		./Common.o\
		./DebugCheck.o\
		./ExonRegions.o\
		./Parsers.o\
		./ReadIndex.o\
		./Sequences.o\
		./SplitAlignment.o\
		./SplitReadAligner.o\
		./ReadStream.o\
		./evalsplitalign.o

splitseq : \
		./AlignmentStream.o\
		./Common.o\
		./DebugCheck.o\
		./ExonRegions.o\
		./ReadIndex.o\
		./Sequences.o\
		./SplitAlignment.o\
		./SplitReadAligner.o\
		./Parsers.o\
		./PathFinder.o\
		./splitseq.o

dosplitaligngen : \
		./AlignmentStream.o\
		./ReadIndex.o\
		./Common.o\
		./DebugCheck.o\
		./ReadIndex.o\
		./Sequences.o\
		./SplitAlignmentGen.o\
		./SplitReadAligner.o\
		./Parsers.o\
		./ReadStream.o\
		./dosplitaligngen.o

decisionmetrics : \
		./Common.o\
		./DebugCheck.o\
		./CompactBreakRegion.o\
		./Parsers.o\
		./decisionmetrics.o

calcspandist : \
		./AlignmentStream.o\
		./Common.o\
		./DebugCheck.o\
		./Parsers.o\
		./calcspandist.o

calccov : \
		./AlignmentStream.o\
		./Common.o\
		./DebugCheck.o\
		./ExonRegions.o\
		./calccov.o

estislands : \
		./EstCatalog.o\
		./estislands.o

localalign : \
		./localalign.o\
		./SimpleAligner.o\
		./DebugCheck.o\
		./Common.o

realigntobreaks : \
		./realigntobreaks.o\
		./Common.o\
		./DebugCheck.o\
		./AlignmentStream.o\
		./Parsers.o\
		./ReadStream.o\
		./Sequences.o\
		./SimpleAligner.o

realign2 : \
		./realign2.o\
		./Common.o\
		./DebugCheck.o\
		./AlignmentProbability.o\
		./AlignmentStream.o\
		./Parsers.o\
		./ReadStream.o\
		./Sequences.o\
		./SimpleAligner.o

realign : \
		./realign.o\
		./Common.o\
		./DebugCheck.o\
		./AlignmentStream.o\
		./Parsers.o\
		./ReadStream.o\
		./Sequences.o\
		./SimpleAligner.o

realign_rna : \
		./realign_rna.o\
		./Common.o\
		./DebugCheck.o\
		./AlignmentProbability.o\
		./AlignmentStream.o\
		./ExonRegions.o\
		./Parsers.o\
		./ReadStream.o\
		./Sequences.o\
		./SimpleAligner.o

realignp : \
		./realignp.o\
		./Common.o\
		./DebugCheck.o\
		./AlignmentProbability.o\
		./AlignmentStream.o\
		./Parsers.o\
		./ReadStream.o\
		./Sequences.o\
		./SimpleAligner.o

testsplit : \
		./Common.o\
		./DebugCheck.o\
		./AlignmentStream.o\
		./ReadStream.o\
		./Sequences.o\
		./testsplit.o

aligntrue : \
		./aligntrue.o\
		./AlignmentStream.o\
		./Common.o\
		./DebugCheck.o\
		./Parsers.o\
		./Sequences.o\
		./SimpleAligner.o

alignnull : \
		./alignnull.o\
		./AlignmentStream.o\
		./Common.o\
		./DebugCheck.o\
		./Parsers.o\
		./ReadStream.o\
		./Sequences.o\
		./SimpleAligner.o

rankclusters : \
		./rankclusters.o\
		./Common.o\
		./DebugCheck.o\
		./Parsers.o

cycles : \
		./Common.o\
		./DebugCheck.o\
		./ShortestPath.o\
		./Parsers.o\
		./cycles.o

paths : \
		./Common.o\
		./DebugCheck.o\
		./ShortestPath.o\
		./Parsers.o\
		./paths.o

filterreads : \
		./AlignmentStream.o\
		./Common.o\
		./DebugCheck.o\
		./RegionDB.o\
		./filterreads.o

filterrnaseq : \
		./AlignmentStream.o\
		./Common.o\
		./DebugCheck.o\
		./ExonRegions.o\
		./filterrnaseq.o

clusterconnectivity : \
		./Common.o\
		./DebugCheck.o\
		./Parsers.o\
		./clusterconnectivity.o

qualtrimfastq : \
		./Common.o\
		./DebugCheck.o\
		./qualtrimfastq.o

testssealign : \
		./Common.o\
		./DebugCheck.o\
		./SimpleAligner.o\
		./testssealign.o

samplefastq : \
		./Common.o\
		./DebugCheck.o\
		./samplefastq.o

bamdiscordantfastq : \
 		./Common.o\
		./DebugCheck.o\
		./bamdiscordantfastq.o\
		$(BAMTOOLS_OBJS)\
		$(ZDIR)/$(ZLIB)

bamsamplespanning : \
		./Common.o\
		./DebugCheck.o\
		./bamsamplespanning.o\
		$(BAMTOOLS_OBJS)\
		$(ZDIR)/$(ZLIB)

bamsamplefastq : \
		./Common.o\
		./DebugCheck.o\
		./bamsamplefastq.o\
		$(BAMTOOLS_OBJS)\
		$(ZDIR)/$(ZLIB)

bamfastq : \
		./Common.o\
		./DebugCheck.o\
		./bamfastq.o\
		$(BAMTOOLS_OBJS)\
		$(ZDIR)/$(ZLIB)

bamconcordantcounts : \
		./Common.o\
		./DebugCheck.o\
		./Parsers.o\
		./bamconcordantcounts.o\
		$(BAMTOOLS_OBJS)\
		$(ZDIR)/$(ZLIB)

bamconcordantreads : \
		./Common.o\
		./DebugCheck.o\
		./Parsers.o\
		./bamconcordantreads.o\
		$(BAMTOOLS_OBJS)\
		$(ZDIR)/$(ZLIB)

bamstats : \
		./Common.o\
		./DebugCheck.o\
		./Parsers.o\
		./bamstats.o\
		$(BAMTOOLS_OBJS)\
		$(ZDIR)/$(ZLIB)

bammixture : \
		./Common.o\
		./DebugCheck.o\
		./bammixture.o\
		$(BAMTOOLS_OBJS)\
		$(ZDIR)/$(ZLIB)
samplegc : \
		./Common.o\
		./DebugCheck.o\
		./Parsers.o\
		./Sequences.o\
		$(BAMTOOLS_OBJS)\
		$(ZDIR)/$(ZLIB)\
		./samplegc.o

samplemap : \
		./Common.o\
		./DebugCheck.o\
		./Parsers.o\
		./Sequences.o\
		$(BAMTOOLS_OBJS)\
		$(ZDIR)/$(ZLIB)\
		./samplemap.o

estimategc : \
		./Common.o\
		./DebugCheck.o\
		./Parsers.o\
		./Sequences.o\
		$(BAMTOOLS_OBJS)\
		$(ZDIR)/$(ZLIB)\
		./estimategc.o

samplemappable : \
		./Common.o\
		./DebugCheck.o\
		./Parsers.o\
		./RegionDB.o\
		./samplemappable.o

clean : 
		rm -f $(BINS) *.o

