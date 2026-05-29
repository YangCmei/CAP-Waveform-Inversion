    implicit none

    integer*4 intorg
    logical*4 logorg
    real*4 realdata,realorg
    character*4 cheader
    common/datr/ realorg(70)
    common/dati/ intorg(35)
    common/datl/ logorg(5)
    common/datd/ realdata(4000001)
    common/datc/ cheader(48)

    character*72 fusename

    fusename='a.sac'
    call readsacb(fusename)

    end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine readsacb(fn)
	character*100 fn
	real*4 realorg
	integer*4 intorg,npts
	logical*4 logorg
	character*4 cheader
	real*4 realdata
	common/datr/ realorg(70)
	common/dati/ intorg(35)
	common/datl/ logorg(5)
	common/datd/ realdata(4000001)
	common/datc/ cheader(48)

	open(10,file=fn,status='old',form='unformatted',recl=1,access='direct')
		do i=1,70
			read(10,rec=i)realorg(i)
		enddo
		do i=1,35
			read(10,rec=i+70)intorg(i)
		enddo
		do i=1,5
			read(10,rec=i+105)logorg(i)
		enddo
		do i=1,48
			read(10,rec=i+110)cheader(i)
		enddo
		npts=intorg(10)
		do i=1,npts
			read(10,rec=i+158)realdata(i)
		enddo
	close(10)
        write(*,*) realorg
        write(*,*) intorg
        write(*,*) logorg
        write(*,*) cheader
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
