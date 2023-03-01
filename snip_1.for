if ( (cnx3 - as) .lt. margin) then	! The East face is pierced for the FIRST TIME
												
    ys = cnz3 + (cnz - cnz3)*(as - cnx3)/(cnx - cnx3)
            

    ! case b.2.1.1  Checking if u(iii,jjj,kkk) is bellow cn3-cn 
    if( ( dom(ib)%zc(kkk) - ys) .lt. margin) then
        if( ( cnx5 - as) .gt. -margin) then	! SAFETy CHECK The EAST face is pierced for the SECOND TIME
            ys2 = cnz5 + (cnz3 - cnz5)*(as - cnx5)/(cnx3 - cnx5)
            if ( ( dom(ib)%zc(kkk) - ys2) .lt. margin) then
                dom(ib)%u(iii,jjj,kkk) = 0.0d00
            endif	
        elseif( ( cnx4 - as) .lt. margin) then	! SAFETy CHECK The EAST face is pierced for by cn4
            ys2 = cnz2 + (cnz4 - cnz2)*(as - cnx2)/(cnx4 - cnx2)
            if ( ( dom(ib)%zc(kkk) - ys2) .lt. -margin) then
                dom(ib)%u(iii,jjj,kkk) = 0.0d00
            endif		

        else
            dom(ib)%u(iii,jjj,kkk) = 0.0d00
            if( ( dom(ib)%zc(kkk)-cnz5) .lt. margin) then
                dom(ib)%u(iii+1,jjj,kkk) = 0.0d00			! There may be an issue
            endif
        endif														
    elseif( ( dom(ib)%zc(kkk)-cnz3) .lt. margin) then	! We do not want go above the top of the geometry
        if( ( dom(ib)%zc(kkk)-cnz5) .lt. margin) then
            dom(ib)%u(iii+1,jjj,kkk) = 0.0d00			! There may be an issue
        endif
    endif



        else 										! cn-cn3 is above u(iii,jjj,kkk)
            dom(ib)%u(iii,jjj,kkk) = 0.0d00
            if ( ( cnx5 - as2 ) .lt. margin ) then
                ys2 = cnz3 + (cnz5-cnz3)*(as2-cnx3)/(cnx5-cnx3)
                if ( ( dom(ib)%zc(kkk-1) - ys2 ) .gt. -margin ) then		
                    dom(ib)%u(iii+1,jjj,kkk) = 0.0d00		
                endif
            else
                dom(ib)%u(iii,jjj,kkk-1) = 0.0d00	
            endif
        endif
    elseif ( ( cnx3 - ae2) .lt. margin) then	! Checking if u(iii+1,jjj,kkk) is above cn2-cn4
        if ( ( cnx6 - ae2 ) .gt. -margin ) then
            ys2 = cnz4 + (cnz6-cnz4)*(ae2-cnx4)/(cnx4-cnx6)
            if ( ( dom(ib)%zc(kkk-1) - ys2 ) .gt. -margin ) then		
                dom(ib)%u(iii+1,jjj,kkk-1) = 0.0d00		
            endif
        else
            dom(ib)%u(iii+1,jjj,kkk-1) = 0.0d00
        endif
    endif
        

    ! case b.2.1.2  Checking if u(iii,jjj,kkk+1) is bellow cn3-cn 
    if( ( dom(ib)%zc(kkk+1) - ys) .lt. margin) then
        if( ( cnx5 - ae) .lt. margin) then	! SAFETy CHECK The EAST face is pierced for the SECOND TIME
            ys2 = cnz5 + (cnz3 - cnz5)*(ae - cnx5)/(cnx3 - cnx5)
            if ( ( dom(ib)%zc(kkk+1) - ys2) .gt. -margin) then
                dom(ib)%u(iii,jjj,kkk+1) = 0.0d00
            endif
        else
            dom(ib)%u(iii,jjj,kkk+1) = 0.0d00
            if( ( dom(ib)%zc(kkk+1)-cnz5) .lt. margin) then
                dom(ib)%u(iii+1,jjj,kkk+1) = 0.0d00			! There may be an issue
            endif
        endif														
    elseif( ( dom(ib)%zc(kkk+1)-cnz3) .lt. margin) then		! We do not want go above the top of the geometry
        if( ( dom(ib)%zc(kkk+1)-cnz5) .lt. margin) then
            dom(ib)%u(iii+1,jjj,kkk+1) = 0.0d00			! There may be an issue
        endif
    endif

    
elseif ( (cnx5 - ae) .gt. -margin) then	! The East face is pierced for the FIRST TIME

    ys = cnz5 + (cnz3-cnz5)*(ae-cnx5)/(cnx3-cnx5)

    ! case b.2.1.0  Checking if u(iii,jjj,kkk+1) is bellow cn3-cn 
    if( ( dom(ib)%zc(kkk+1) - ys) .lt. margin) then
        if( ( dom(ib)%zc(kkk+1)-cnz5) .lt. margin ) then
            dom(ib)%u(iii,jjj,kkk+1) = 0.0d00
        endif														
    endif


    ! case b.2.1.1  Checking if u(iii,jjj,kkk) is bellow cn3-cn 
    if( ( dom(ib)%zc(kkk) - ys) .lt. margin) then
        if( ( cnx4 - ae) .gt. -margin) then	! SAFETy CHECK The EAST face is pierced for by cn4
            ys2 = cnz2 + (cnz4 - cnz2)*(ae - cnx2)/(cnx4 - cnx2)
            if ( ( dom(ib)%zc(kkk) - ys2) .gt. -margin) then
                dom(ib)%u(iii,jjj,kkk) = 0.0d00
            endif														
        elseif( ( dom(ib)%zc(kkk)-cnz5) .lt. margin ) then
            dom(ib)%u(iii,jjj,kkk) = 0.0d00														
        endif
    endif															


else

    flag2 = 1

endif !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	Finished Check of East Face


if(flag1.eq.1) then
    if(flag2.eq.1) then

        if ( ( dom(ib)%zc(kkk-1) - cnz4) .gt. -margin) then
            dom(ib)%u(iii,jjj,kkk-1) = 0.0d00													
        endif

        dom(ib)%u(iii,jjj,kkk) = 0.0d00

        if ( ( dom(ib)%zc(kkk+1) - cnz3) .lt. margin) then
            dom(ib)%u(iii,jjj,kkk+1) = 0.0d00																									
        endif

    endif
endif !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	Finished Additional Check of East Face
