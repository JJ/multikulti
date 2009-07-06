# Multikulti algorithm: an island model implementation of a GA>
#     Copyright (C) <2009>  <Lourdes Araujo, Juan Julian Merelo>

#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>

PVM_ROOT	=	/usr/share/pvm3
SDIR	=	.
#BDIR	=	$(HOME)/pvm3/bin
#BDIR	=	$(SDIR)
#BDIR	=	/usr/bin/
BDIR	=	/usr/share/pvm3
XDIR	=	$(SDIR)
ARCHLIB	=	

CC	=	c++
OPTIONS	=	-O
CFLAGS	=	$(OPTIONS) -I$(PVM_ROOT)/include $(ARCHCFLAGS)
LIBS	=	-lpvm3 $(ARCHLIB)
GLIBS	=	-lgpvm3


LFLAGS	=	-L$(PVM_ROOT)/lib/$(PVM_ARCH)

default:	control funciones

clean:
	rm -f *.o control funciones 

$(XDIR):
	- mkdir $(BDIR)
	- mkdir $(XDIR)

control:	control.o 
	$(CC) $(CFLAGS) -Wno-deprecated -o control control.o $(LFLAGS) $(LIBS)


funciones:	funciones.o
	$(CC) $(CFLAGS) -Wno-deprecated -o funciones funciones.o $(LFLAGS) $(LIBS) 

control.o:	control.c	        
	$(CC) -Wno-deprecated -c control.c 

funciones.o:	funciones.c	        
	$(CC) -Wno-deprecated -c funciones.c 

