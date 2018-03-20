	INTEGER FUNCTION NARGU(ARG)

c	Count how many time steps have been introduced in string arg

c	Variables
	character*(*) arg,buf*1024

c	Procedure
	i = 0
	l=largo(arg)
	if (l.eq.0) then
		nargu=0
		return
	endif
	if (arg(l:l).eq.',') arg(l:l)=' '
	buf=arg
	do while (index(buf,',').gt.0)
	i=i+1
	buf=buf(index(buf,',')+1:)
	enddo
	nargu=i+1
	return
	end
