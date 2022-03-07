
CC = cc
#FLAGS = -g
OPENWINHOME = /usr/openwin
MFILE = 
INCDIR = -I$(OPENWINHOME)/include
LIBDIR = -L$(OPENWINHOME)/lib
LIBS = -lxview -lolgx -lX11 

libs.o = Alloc.o HGLfuncs.o

all: mapview MakeCons Consto01mask PrintStrat Translate heapsortHGL DotPlotTool

Alloc.o: Alloc.c 
	$(CC) $(FLAGS) -c Alloc.c

HGLfuncs.o: HGLfuncs.c 
	$(CC) $(FLAGS) -c HGLfuncs.c

MAP_ChooseFile.o: MAP_ChooseFile.c 
	$(CC) $(FLAGS) -c MAP_ChooseFile.c $(INCDIR)

mapview: mapview.c MAP_ChooseFile.o $(libs.o) 
	$(CC) $(FLAGS) -o $@ $@.c MAP_ChooseFile.o $(libs.o) $(INCDIR) $(LIBDIR) $(LIBS)

DotPlotTool: DotPlotTool.c MAP_ChooseFile.o $(libs.o) 
	$(CC) $(FLAGS) -o $@ $@.c MAP_ChooseFile.o $(libs.o) \
$(INCDIR) $(LIBDIR) $(LIBS)

MakeCons: MakeCons.c $(libs.o) $(MKFILE)
	$(CC) $(FLAGS) -o $@ $@.c $(libs.o)

Translate: Translate.c $(libs.o) $(MKFILE)
	$(CC) $(FLAGS) -o $@ $@.c $(libs.o)

heapsortHGL: heapsortHGL.c $(libs.o) $(MKFILE)
	$(CC) $(FLAGS) -o $@ $@.c $(libs.o)

PrintStrat: PrintStrat.c $(libs.o) $(MKFILE)
	$(CC) $(FLAGS) -o $@ $@.c $(libs.o)

Consto01mask: Consto01mask.c $(libs.o) $(MKFILE)
	$(CC) $(FLAGS) -o $@ $@.c $(libs.o)
