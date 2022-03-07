OBJS= ParseMenu.o main.o BasicDisplay.o EventHandler.o FileIO.o \
DrawNA.o Free.o BuiltIn.o Edit.o Genbank.o Scroll.o ChooseFile.o \
CutCopyPaste.o HGLfile.o
SRCS= ParseMenu.c main.c BasicDisplay.c EventHandler.c FileIO.c \
DrawNA.c Free.c BuiltIn.c Edit.c Genbank.c Scroll.c ChooseFile.c \
CutCopyPaste.c HGLfile.c

LIBS= -lm -lxview -lolgx -lX11
CFLAGS= -g -L/usr/openwin/lib -I/usr/openwin/include
CC = cc
#	Possible defines, SUN4 SGI DEC HGL
DEFINES =  -DLINUX

GLOBAL_DEPENDS= defines.h menudefs.h

gde: $(OBJS)
	$(CC) -o $@ $(OBJS)  $(LIBS)  $(CFLAGS) $(DEFINES)

ParseMenu.o: ParseMenu.c $(GLOBAL_DEPENDS)
	$(CC) -c ParseMenu.c   $(CFLAGS) $(DEFINES)

main.o: main.c $(GLOBAL_DEPENDS) globals.h 
	$(CC) -c main.c  $(CFLAGS) $(DEFINES)

BasicDisplay.o: BasicDisplay.c $(GLOBAL_DEPENDS)
	$(CC) -c BasicDisplay.c  $(CFLAGS) $(DEFINES)

EventHandler.o: EventHandler.c $(GLOBAL_DEPENDS)
	$(CC) -c EventHandler.c  $(CFLAGS) $(DEFINES)

FileIO.o: FileIO.c $(GLOBAL_DEPENDS)
	$(CC) -c FileIO.c $(CFLAGS) $(DEFINES)

DrawNA.o: DrawNA.c $(GLOBAL_DEPENDS)
	$(CC) -c DrawNA.c $(CFLAGS) $(DEFINES)

Edit.o: Edit.c $(GLOBAL_DEPENDS)
	$(CC) -c Edit.c $(CFLAGS) $(DEFINES)

Free.o: Free.c $(GLOBAL_DEPENDS)
	$(CC) -c Free.c $(CFLAGS) $(DEFINES)

BuiltIn.o: BuiltIn.c $(GLOBAL_DEPENDS)
	$(CC) -c BuiltIn.c $(CFLAGS) $(DEFINES)

Genbank.o: Genbank.c $(GLOBAL_DEPENDS)
	$(CC) -c Genbank.c $(CFLAGS) $(DEFINES)

Scroll.o: Scroll.c $(GLOBAL_DEPENDS)
	$(CC) -c Scroll.c $(CFLAGS) $(DEFINES)

ChooseFile.o: ChooseFile.c $(GLOBAL_DEPENDS)
	$(CC) -c ChooseFile.c $(CFLAGS) $(DEFINES)

CutCopyPaste.o: CutCopyPaste.c $(GLOBAL_DEPENDS)
	$(CC) -c CutCopyPaste.c $(CFLAGS) $(DEFINES)

HGLfile.o: HGLfile.c $(GLOBAL_DEPENDS)
	$(CC) -c HGLfile.c $(CFLAGS) $(DEFINES)

