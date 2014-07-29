        program monthly alpha
        implicit none

*************************************************************
* A programme to predict monthly alpha based on DavidsPAR
* program. Therefore the main program is the same as DavidsPAR
* file, except for monthly alpha calculation and output are
* added.
* Repeating might be involved if the soil moisture was not
* equal at the end of year
* negative alpha at high latitude in winter time is assigned as 
* zero value.
*********************************************************
* The input data is monthly temperature, precipitation and
* cloudness. Soil field capacity is needed as well, and a 
* global value of 150 mm is usually be used
*********************************************************
*Variables used:
* maxgrd = maximum number of  cells, set to 259200, because 
* lon:720 grids, lat=360 grids (-90 to 90 degrees) 
* sysout = standard output device (console, screen)
* sysin  = standard input device  (console, keyboard)
*********************************************************

        integer maxgrd,sysout,sysin
        parameter(maxgrd=571290)
* input and output devices for vax/vms and ibm/pc-at
        parameter(sysout=6,sysin=5)
      

* common for environmental variables and multipliers.
   
        real    gdd(2006),gsdri(maxgrd),m4dri(maxgrd),tcold(maxgrd)
        real    twarm(maxgrd),tmannual(maxgrd),childy(maxgrd)
        real    nodata,yeardri(maxgrd),yeararid(maxgrd),gdd0(maxgrd)
        real    gdd05(maxgrd),gdd10(maxgrd)
        real    par0(maxgrd),par5(maxgrd),yrpar(maxgrd)

        real	yralpha(maxgrd),mpar0(maxgrd),mgdd0(maxgrd)
        real    mgdd5(maxgrd),mgdd10(maxgrd)

        parameter(nodata=-99.8)

        common  /climo/gdd,gdd0,gdd05,gdd10,gsdri,m4dri,yeardri,
     %          yeararid,tcold,twarm,tmannual,childy
    
* local variables

        integer i,k,ll,ncell
        real    lat,lon,plat(maxgrd),plon(maxgrd),clon(maxgrd)
        real    clat(maxgrd),mtc(12),mpr(12),mcl(12),fcap
        real    imcl(maxgrd,12),impr(maxgrd,12),inputmtc(12)
        real    minlat,maxlat,minlon,maxlon
        character*1 ans
        logical blank

* open output file

        open(unit=3,file='out_mPAR_0C.txt',status='unknown')
        open(unit=4,file='out_PAR_0C.txt',status='unknown')
        open(unit=15,file='out_GDD_0C.txt',status='unknown')
        open(unit=16,file='out_mGDDay_0C.txt',status='unknown')
        open(unit=17,file='out_MI.txt',status='unknown')
        open(unit=18,file='out_yralpha.txt',status='unknown')
        open(unit=19,file='out_malpha.txt',status='unknown')

        open(unit=30,file='out_GDD_5C.txt',status='unknown')
        open(unit=21,file='out_mGDDay_5C.txt',status='unknown')
        open(unit=22,file='out_GDD_10C.txt',status='unknown')
        open(unit=23,file='out_mGDDay_10C.txt',status='unknown')
* Get field capacity
       
        write(sysout,'(/''  Give field capacity (mm): ''$)')
        read(sysin,*)fcap

* call and open climate input file (this needs redoing)

        open(10,file='in_temp.txt',status='old')
        open(11,file='in_rain.txt',status='old')
        open(12,file='in_sun.txt',status='old')


* Check if full boreal  window needed or some other size

        write(sysout,'(/'' Use Australian coords? (y/n): ''$)')
        read(sysin,'(a1)')ans
        if(ans.eq.'n'.or.ans.eq.'N')then
            write(sysout,'(/'' Use global lat,-90 to 90 deg''$)')
            write(sysout,'(/'' Give minimum latitude: ''$)')
            read(sysin,*)minlat
            write(sysout,'(/'' Give maximum latitude: ''$)')
            read(sysin,*)maxlat
            write(sysout,'(/'' Give minimum longitude: ''$)')
            read(sysin,*)minlon 
            write(sysout,'(/'' Give maximum longitude: ''$)')
            read(sysin,*)maxlon
        else
*          minlat=-90.
*          maxlat=90.
*          minlon=-180.
*          maxlon=180.
            minlat=-43.725
            maxlat=-9.025
            minlon=112.925
            maxlon=153.975
        endif

        write(sysout,'(/'' Selected window has the co-ordinates:'')')
        write(sysout,'(/'' minimum latitude= '',f8.3,
     %                  /'' maximum latitude= '',f8.3,
     %                  /'' minimum longitude='',f8.3,
     %                  /'' maximum longitude='',f8.3)')minlat,maxlat,
     %                                                minlon,maxlon
*Initialize parameters and start loop to read and call subroutines
        ncell=0    
100     do 200 ll=1,maxgrd
* Read in precip and cloudiness data into an array and then check for lon and 
* lat with temperature file
            read(11,*,err=999,end=900)plon(ll),plat(ll),
     %          (impr(ll,i),i=1,12)
            read(12,*,err=999,end=900)clon(ll),clat(ll),
     %          (imcl(ll,i),i=1,12)
200     continue

        write(sysout,'(/'' Precip. data read into an array'')')
        write(sysout,'(/'' Cloudiness data read into an array'')')

        write(sysout,'(/'' Calculate for all cells '')')
* Give each gridcell a number and call main subroutine for each cell. Main cycle 
        do 201 ll=1,maxgrd
        print *, 'Completion: ', real(ll)/real(maxgrd)*100, '% Done'
* Checking on lats and longs of each cell and calling main routine if OK
            read(10,*,err=999,end=900)lon,lat,(inputmtc(i),i=1,12)
            do 202 i=1,12
                mtc(i)=inputmtc(i)
202         continue
* check first to see if within window
        if(lon.ge.minlon.and.lon.le.maxlon
     %      .and.lat.ge.minlat.and.lat.le.maxlat)then
* find correct precip cell
            do 205 k=1,maxgrd
                if(lon.eq.plon(k).and.lat.eq.plat(k))then 
                    do 206 i=1,12
                    mpr(i)=impr(k,i)
206                 continue
                goto 210
                endif
205         continue
* Find correct cloudiness cell
210         do 215 k=1,maxgrd
                if(lon.eq.clon(k).and.lat.eq.clat(k))then
                    do 207 i=1,12
                        mcl(i)=imcl(k,i)
207                 continue
                    goto 220
                endif
215         continue
220         ncell=ncell+1        
* call environmental subroutines, count cells, flag lat & lon to be used
            if(mpr(1).gt.nodata.or.mtc(1).gt.nodata)then
                call staenv(sysout,lon,lat,mtc,mpr,mcl,ll,blank,fcap)
            else
                write(19,'(2f9.2,12f14.1)')lat,lon,nodata,nodata,
     %		        nodata,nodata,nodata,nodata,nodata,nodata,
     %		        nodata,nodata,nodata,nodata
            endif
        endif               
201     continue

900     write(sysout,'(i6,'' cells processed'')')ncell                  
        stop 'KLIMAT finished'

999     write(sysout,300)
        stop 'PEATKLIMAT terminated because of read errors.'     
2       format(x,2(f8.2),f8.2)
*90     format(e18.7,e18.7,12e18.7)
*91     format(e18.7,e18.7,12e18.7)
*92     format(e18.7,e18.7,12e18.7)
300     format(/' error in reading environmental data')

        end program

******************************************************************************* 

        subroutine staenv(sysout,lon,lat,mtc,mpr,mcl,ll,blank,fcap)
        implicit none 

* gets and prepares environmental data for use in peatstash
* gets climate from gridded datasets for temp, precip, sunshine 
* works out gdd, drought index,mean annual temperature 
* light intensity, temp coldest month, temp. warmest month
* written by martin sykes january 1992, adapted Angela09 for peatl 

       integer sysout,ll
       integer sysin,maxgrd
       real lon,lat,mtc(12),mpr(12),mcl(12)
       logical blank

       parameter (sysin=5,maxgrd=571290)

* common environmental variables
* 
* gdd    - growing degree days - site
* childy - number of chill days with temp <5 degrees
* gsdri  - growing season drought index - site
* m4dri  - -4 oc drought index - site
* yeardri - all year round drought index- site
* yeararid - all year round aridity index- site
* tcold  - temperature coldest month  - site
* twarm  - temperature warmest month  - site
* tmmanual - mean annual temperature
 
        real	mgdd0(maxgrd),mgdd5(maxgrd),mgdd10(maxgrd)
        real    mpar0(maxgrd),yralpha(maxgrd)
        real    gdd(maxgrd),gsdri(maxgrd),m4dri(maxgrd),tcold(maxgrd)
        real    twarm(maxgrd),childy(maxgrd),tmannual(maxgrd)
        real    gdd0(maxgrd),gdd05(maxgrd),gdd10(maxgrd)
        real    yeararid(maxgrd),par5(maxgrd),par0(maxgrd),yrpar(maxgrd)
        real    ttot(maxgrd),gprec(maxgrd),gevap(maxgrd),yeardri(maxgrd)
        integer gdday,gdday5,gdday10
        common  /climo/gdd,gdd0,gdd05,gdd10,gsdri,m4dri,yeardri,
     %          yeararid,Tcold,twarm,tmannual,childy
* local variables
        integer k
        real    fcap
* set blank flag to false before calling env routines
        blank=.false.
        call stasub(sysout,lat,lon,mtc,mpr,mcl,ll,blank,fcap,
     %              gprec,gevap,par0,par5,yrpar,
     %              gdday,gdday5,gdday10,
     %              mgdd0,mgdd5,mgdd10,
     %              yralpha,mpar0)
* if this cell is a cell with no information possible go back to main routine
        if(blank) write(sysout,'(/'' Cell with no information'')')
        if(blank) return
* else work out which is temperature of coldest month
* and warmest month for cell
        tcold(ll)=50.0
        twarm(ll)=-50.0
        ttot(ll)=0.0
        do 20 k=1,12
            if(mtc(k).lt.tcold(ll)) tcold(ll)=mtc(k)
            if(mtc(k).gt.twarm(ll)) twarm(ll)=mtc(k)
20      continue

        do 30 k=1,12
            ttot(ll)=ttot(ll)+mtc(k)
30      continue
        tmannual(ll)=ttot(ll)/12
        write(15,102)lat,lon,gdd0(ll)
        write(3,102)lat,lon,mpar0(ll)
        write(4,102)lat,lon,par0(ll)
        write(16,*)lat,lon,mgdd0(ll)
        write(18,102)lat,lon,yralpha(ll)
        write(17,102)lat,lon,yeararid(ll)
        write(30,102)lat,lon,gdd05(ll)
        write(21,*)lat,lon,mgdd5(ll)
        write(22,102)lat,lon,gdd10(ll)
        write(23,*)lat,lon,mgdd10(ll)
* Format statements
102     format(2f9.2,f14.6)
        
        return

	end subroutine
     
*******************************************************************************

        subroutine stasub(sysout,lat,lon,mtc,mpr,mcl,ll,blank,fcap,
     %                  gprec,gevap,par0,par5,yrpar,
     %                  gdday,gdday5,gdday10,
     %					mgdd0,mgdd5,mgdd10,yralpha,mpar0)
        implicit none 
        integer maxgrd,year,sysout
        real tref
        parameter (maxgrd=571290,year=365,tref=5.0)

* originally written by wolfgang cramer
* substantially amended for use in forska2 - climate and landscape versions 
* by colin prentice and martin sykes 1990
* further amended to include new evapo routine september 1991 martin sykes
* later temporary fix for insolation, minor health warning december 1991
* martin sykes, uppsala.
* this version for use in stash january 1992 mts
* then adapted for use in bioclimatic envelope modeling of peatlands
*
* common environmental variables 
*
* gdd    - growing degress days  - site 
* childy _ number of chill days c < 5 degrees c - site
* gsdri  - growing season drought index - site
* m4dri  - -4 oc drought index - site
* yeardri - all year round drought index - site
* gsins  - growing season average light intensity - site (dpar)
* m4ins  - -4 oc season average light intensity - site (dpar)

        real    mgdd0(maxgrd),mgdd5(maxgrd),mgdd10(maxgrd)
        real    yralpha(maxgrd),mpar0(maxgrd)
        real    gdd(maxgrd),gsdri(maxgrd),m4dri(maxgrd),tcold(maxgrd)
        real    twarm(maxgrd),childy(maxgrd),gprec(maxgrd),gevap(maxgrd)
        real    yeardri(maxgrd),yeararid(maxgrd)
        real    par5(maxgrd),par0(maxgrd)
        real    yrpar(maxgrd),gdd0(maxgrd),gdd05(maxgrd),gdd10(maxgrd)
        common  /climo/gdd,gdd0,gdd05,gdd10,gsdri,m4dri,yeardri,
     %          yeararid,tcold,twarm,childy
* local variables
        integer days(12),dn,ind,j,k,ll,runc,sta,int,m4day,gdday,yday
        integer gdday5,gdday10

        real    clou(365),dsm(365),prec(365),temp(365)
        real	mcl(12),mpet(12),maet(12),malpha(12),mpr(12),mtc(12),dpar
        real	alpha,cw,dpet,eccen,k1,k2,k3,lat,lon,lsm
        real	solc,spl,ysm,daet,rno,foudpt,foudae,tgsdpt,tgsdae
        real    fcap,yraet,yrpet,meanprec,yearprec,dpar0,dpar5,dyrpar
        integer out
        logical lr,blank
        parameter(out=11)
        parameter(ind=10,sta=14,alpha=0.17,solc=1360.,cw=1.,
     %		    k1=610.78,k2=17.269,k3=237.3,eccen=0.01675)
        data    days/31,28,31,30,31,30,31,31,30,31,30,31/
* make daily arrays of things we have 
        call daily(mtc,temp) 
        call daily(mpr,prec) 
        call daily(mcl,clou)
        dn=0
        yearprec=0.
        yrpet=0.
        yraet=0.
        do 500 j=1,12
            do 500 k=1,days(j)
                dn=dn+1
                prec(dn)=prec(dn)/days(j)
500         continue	
* daynumber and yearly totals are set to zero
        dsm(1)=fcap
        lsm=fcap
        lr=.true.
        runc=1
* jump back here on rerun (if soil moisture at end of year was not equal
* the value it began with)
511     dn=0
        rno=0.
        foudpt=0.
        foudae=0.
        tgsdpt=0.
        tgsdae=0.
        dpar5=0.0
        dpar0=0.0
        dyrpar=0.0
        m4day=0
        gdday=0
        gdday5=0
        gdday10=0
        gdd(ll)=0.0
        gdd0(ll)=0.0
        gdd05(ll)=0.0
        gdd10(ll)=0.0
        childy(ll)=0.0
        mgdd0(ll)=0.0
        mgdd5(ll)=0.0
        mgdd10(ll)=0.0
        mpar0(ll)=0.0
        yralpha(ll)=0.0
*
* ------------------ start monthly calculation loop --------------------
*
        do 512 j=1,12
            mpet(j)=0.
            maet(j)=0.
*
* ------------------- start daily calculation loop --------------------
*
            do 513 k=1,days(j)
                dn=dn+1
                if(dn.eq.1) then
                    ysm=lsm
                else
                    ysm=dsm(dn-1)
                endif
                spl=cw*ysm/fcap
                call evapo(alpha,eccen,solc,k1,k2,k3,lat,temp(dn),
     %   	             clou(dn),dn,spl,daet,dpet,dpar)
*
* the soil receives today's amount of rain and it loses today's amount of
* "actual" evapotranspiration
                dsm(dn)=ysm+prec(dn)-daet
* adjust for bucket size
                if(dsm(dn).gt.fcap) then
                    rno=rno+(dsm(dn)-fcap)
                    dsm(dn)=fcap
                else
                    if(dsm(dn).lt.0) dsm(dn)=0.0
                endif
* sum up for monthly or annual values
                mpet(j)=mpet(j)+dpet
                maet(j)=maet(j)+daet
* is this a growing season day or a -4 day ? 
* if so then add to growing days and excess temp to gdd total for year
* also total dpet,daet and dpar for both growing season and -4 season
* total daily light intensity over each season
                if(temp(dn).ge.0.0)then
                    m4day=m4day+1
                    foudpt=foudpt+dpet
                    foudae=foudae+daet
                    dpar0=dpar0+dpar
                    gdday=gdday+1
                    gdd0(ll)=gdd0(ll)+(temp(dn)-0.0)
                endif
                if( temp(dn).ge.5.0) then
                    gdday5=gdday5+1
                    gdd05(ll)=gdd05(ll)+(temp(dn)-5.0)
                endif
                if( temp(dn).ge.10.0) then
                    gdday10=gdday10+1
                    gdd10(ll)=gdd10(ll)+(temp(dn)-10.0)
                endif
                if(temp(dn).gt.-400)then
                    yday=yday+1
                    yrpet=yrpet+dpet
                    yraet=yraet+daet
                    dyrpar=dyrpar+dpar
                endif
*                          if(k.eq.21.and.j.eq.6)then
*                             write(16,105)lon,lat,dpar
*                          end if
*
* here we count the number of chill days i.e c - the number of days with
* mean temp < 5
                if(temp(dn).lt.5.0) childy(ll)=childy(ll)+1
513         continue

*
* -------------------- end of daily loop ------------------------------
*
            malpha(j)=maet(j)/mpet(j)
            if(malpha(j).lt.0.)then
                malpha(j)=0.
            endif
			
512     continue
*
* ------------------------ end of monthly loop ----------------------------

* work out drought index for both growing season and -4 season
* average daily light intensity for both seasons
* check if this cell has anything that can be plotted if not go back

        if(yrpet.le.0)then
          blank=.true.
          return
        endif

*Aridity calculation needs yearly precipitation: 

        do 450 j=1,12
      	  yearprec=yearprec+mpr(j)
450     continue
              

        gsdri(ll)=(tgsdpt-tgsdae)/tgsdpt
        m4dri(ll)=(foudpt-foudae)/foudpt
        yeardri(ll)=(yrpet-yraet)/yrpet
        yeararid(ll)=yearprec/yrpet
        gprec(ll)=tgsdpt
        gevap(ll)=tgsdae
        par0(ll)=dpar0/1000000.0
        par5(ll)=dpar5/1000000.0
        yrpar(ll)=dyrpar/1000000.0
        yralpha(ll)=yraet/yrpet
        mpar0(ll)=par0(ll)/gdday
        mgdd0(ll)=gdd0(ll)/gdday
        mgdd5(ll)=gdd05(ll)/gdday5
        mgdd10(ll)=gdd10(ll)/gdday10

        if(ll.lt.5) print*,lon,lat,yrpet,yraet,yearprec

*        write(17,*)lon,lat,yearprec,yrpet,yeararid(ll)
* test for soil moisture replenishment

        if(int(lsm*10).ne.int(dsm(365)*10)) then
* changed to 30 sept 1992 
            if(runc.gt.40) then
                write(*,1201) dsm(365)
1201            format(' stability still not reached..:',f8.3)
            endif
            lsm=dsm(365)
            lr=.false.
            runc=runc+1
105         format(f6.1,f6.1,f14.6)
            goto 511
        endif
        continue
        write(19,103)lat,lon,(malpha(j),j=1,12)
103     format(2f9.2,12f14.6)

        end subroutine

*---------------------------------------------------------------------------*

        subroutine daily(mly,dly)

        implicit none

        integer days(14),daft,dbef,dn,j,k
        real    dly(365),mly(12),tmly(14),d15,inc,maft,mbef
        data    days/31,31,28,31,30,31,30,31,31,30,31,30,31,31/
                Tmly(1)=mly(12)

        do 490 j=1,12
            tmly(j+1)=mly(j)
490     continue
        tmly(14)=mly(1)
        dn=0
        do 500 j=2,13
            do 500 k=1,days(j)
                dn=dn+1
*  first half of month
                if(k.lt.(days(j)/2.)) then
                    mbef=tmly(j-1)
                    maft=tmly(j)
                    dbef=days(j-1)
                    daft=days(j)
                else
*  second half of month
                    mbef=tmly(j)
                    maft=tmly(j+1)
                    dbef=days(j)
                    daft=days(j+1)
                endif
                inc=(maft-mbef)/((dbef+daft)/2.)
                if(k.lt.(days(j)/2.)) then
                    d15=k+(dbef/2)
                else
                    d15=k-(daft/2)
                endif
                dly(dn)=mbef+(d15*inc)
500     continue

        end subroutine

*-----------------------------------------------------------------------------*
        subroutine evapo(alpha,eccen,solc,k1,k2,k3,lat,dtc,
     %                      clou,day,spl,daet,dpet,dpar)

* this version added into forska2 to correct faults, but dpar calculation
* added from earlier version september 1991. martin sykes
*
* calculates evaporation
* from delta (solar declination, based on day (daynumber)), lat (latitude),
* clou/dcl (cloudiness), sat (slope of water vaporization curve), gamma
* (psychrometer constant), lambda (latent heat of water vaporization),
* solc (solar constant), alpha (albedo), dtc (temperature), spl (supply
* function)
*
* nb: requires now cloudiness as ratio (<1.0)! (i.e dcl=clou/100 (mts))
* Added 2009AGS: d2r= to convert degrees into radians (because gfortran does not accept 
* cosd, sind or tand). 
*
        implicit none

        integer day
        real    alpha,arg,daet,clou,dcl,dcon,ddtc,delta,dpet,dtc,eccen
        real    dpar,gamma,h0,h1,k1,k2,k3,lambda,lat,msolc,ndtc,d2r
        real    hs,pi,sat,solc,spl,u,v,x,y

        pi=4*atan(1.)
        d2r=pi/180.
        dcl=clou/100
        arg=360.*day/365.
        msolc=solc*(1.+2.*eccen*cos(arg*d2r))
* delta is today's solar declination in degrees
        delta=-23.4*cos((((day+10.)*360.)/365.)*d2r)
* gamma, the psychrometer constant (pa/oc), and lambda, the latent heat of
* water vaporisation (mj/kg) are taken from a table
        call table(dtc,gamma,lambda)
* find for the present temperature: sat (pa/oc)
        sat=k1*k2*k3*exp(k2*dtc/(k3+dtc))/((k3+dtc)**2)
* h0 is the point in time when net radiation crosses zero
        u=(0.0036/lambda)*(sat/(sat+gamma))*
     %      (msolc*(0.25+0.5*dcl)*(1.-alpha)*sin(lat*d2r)*sin(delta*d2r)
     %      -((0.2+0.8*dcl)*(107-dtc)))
        v=(0.0036/lambda)*(sat/(sat+gamma))*
     %      (msolc*(0.25+0.5*dcl)*(1-alpha)*cos(lat*d2r)*cos(delta*d2r))
        if((spl-u).ge.v) then
*  supply exceeds demand
            h1=0.
            elseif((spl-u).le.(0.-v)) then
*  demand exceeds supply
            h1=pi
            h0=pi
        else
            arg=((spl-u)/v)
            h1=acos(arg)
        endif
        if(u.ge.v) then
*  polar day
            h0=pi
            elseif(u.le.(0.-v)) then
*  polar night
            h1=0.
            h0=0.
        else
*  normal day and night
            arg=(u/v)*(-1.)
            h0=acos(arg)
        endif
        ddtc=dtc
        ndtc=dtc
        u=(0.0036/lambda)*(sat/(sat+gamma))*
     %      (msolc*(0.25+0.5*dcl)*(1.-alpha)*sin(lat*d2r)*sin(delta*d2r)
     %      -((0.2+0.8*dcl)*(107.-ddtc)))
        dcon=(0.0036/lambda)*(sat/(sat+gamma))*(0.2+0.8*dcl)*
     %      (107.-ndtc)
        daet=((2.*u*(h0-h1)+2.*v*(sin(h0)-sin(h1))+2.*spl*h1)/(pi/12.))-
     %      dcon
        dpet=((2.*u*h0)+(2.*v*sin(h0)))/(pi/12.)
* this is the dpar bit taken from old evapo routines corrected for lat
*  - lambda cock-up september 1991.
* old clpr was clou/100 which is new dcl
        x=(0.25+0.5*dcl)*(1.-alpha)*sin(lat*d2r)*sin(delta*d2r)
        y=(0.25+0.5*dcl)*(1.-alpha)*cos(lat*d2r)*cos(delta*d2r)
        if(x.gt.y) then
            hs=pi
            elseif(x.lt.(-1.*y)) then
            hs=0.
        else
            hs=acos(-1.*tan(lat*d2r)*tan(delta*d2r))
        endif
* calculate total daily par required for sum over the growing season 
        dpar=(solc/0.22)*(0.25+0.5*dcl)*(1.-alpha)*86400*
     %      ((hs*sin(lat*d2r)*sin(delta*d2r))+
     %  	(cos(lat*d2r)*cos(delta*d2r)*sin(hs)))/
     %  	(2*pi)

        end subroutine

*---------------------------------------------------------------------------*

        subroutine table(tc,gamma,lambda)

        implicit none

        integer i,k,l
        real    gbase(2,11),lbase(2,11),gamma,lambda,tc
        data ((gbase(k,l),k=1,2),l=1,11)/-5.,64.6,0.,64.9,5.,65.2,10.,
     %	     65.6,15.,65.9,20.,66.1,25.,66.5,30.,66.8,35.,67.2,40.,67.5,
     %       45.,67.8/
        data ((lbase(k,l),k=1,2),l=1,11)/-5.,2.513,0.,2.501,5.,2.489,
     %  	10.,2.477,15.,2.465,20.,2.454,25.,2.442,30.,2.43,35.,2.418,
     %  	40.,2.406,45.,2.394/

        do 500 i=1,11
* temperature at or below loop value - set loop gamma and return
            if(tc.le.gbase(1,i)) then
                gamma=gbase(2,i)
                lambda=lbase(2,i)
                return
* temperature above highest loop value - set highest gamma and return
            elseif(tc.gt.gbase(1,11)) then
                gamma=gbase(2,11)
                lambda=lbase(2,11)
                return
            endif
500     continue

        end subroutine

*---------------------------------------------------------------------------*
