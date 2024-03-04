      SUBROUTINE shine(tfsn, tfsg, ts, hgbq, hnbq, alb_c, alb_o)
!-----------------------------------------------------------------------------!
!
!  Computes albedo of snow-sea ice following SHINE &
!  HENDERSSON-SELLERS [1985].
!  Simplified by Martouf for use in LIM1D
!  fucking scientists!!!
!
      INCLUDE 'type.com'
!     INCLUDE 'para.com'
!     INCLUDE 'const.com'
!     INCLUDE 'bloc.com'
!     INCLUDE 'ice.com'
!     INCLUDE 'forcing.com'
 
!  OUTPUTS :
!  alb_c : Albedo for clear sky.
!  alb_o :  Albedo for overcast sky.
!
!-----------------------------------------------------------------------------!
!  1) Computation of surface albedo.                                          !
!-----------------------------------------------------------------------------!
!
      WRITE(84,*) ' shine '
      WRITE(84,*) ' ~~~~~~'

      cgren  = 0.06 ! corr factor under cloudy skies (Grenfell Perovich 84)
      alphd  = 0.80 ! fixed boundary values of albedo
      alphdi = 0.72
      alphs  = 0.65
      albice = 0.53 ! albedo of melting ice
!
      if (hnbq.gt.0.0) then                                       
      !------------------
      ! Snow covered-ice
      !------------------
        if (ts.lt.tfsn) then                                      

      ! Cold snow.
      !------------
          if (hnbq.gt.0.05) then                                 
            alb_c = alphd                                        
          else                                                      
            if (hgbq.gt.1.5) then                                
              alb_c = alphdi+(hnbq*(alphd-alphdi)/0.05)          
            else if (hgbq.gt.1.0.and.hgbq.le.1.5) then         
                   al = 0.472+2.0*(alphdi-0.472)*(hgbq-1.0)
            else if (hgbq.gt.0.05.and.hgbq.le.1.0) then        
                   al = 0.2467+(0.7049*hgbq)-(0.8608*(hgbq*hgbq))+
     &                 (0.3812*(hgbq*hgbq*hgbq))                     
            else                                                    
     &        al = 0.1+3.6*hgbq                                  
            endif                                                   
            if (hgbq.le.1.5) alb_c=al+(hnbq*(alphd-al)/0.05)
          endif                                                     
        else                                                        
      ! Melting snow.
      !--------------
          if (hnbq.ge.0.1) then                                  
            alb_c = 0.65                                           
            alb_c = alphs                                      
          else                                                      
            alb_c = albice+((alphs-albice)/0.1)*hnbq
          endif                                                     
        endif                                                       
      else                                                          
      !----------
      ! Bare ice
      !----------
        if (ts.lt.tfsg) then                                      

      ! Cold ice
      !-----------
          if (hgbq.gt.1.5) then                                  
            alb_c = alphdi                                          
          else if (hgbq.gt.1..and.hgbq.le.1.5) then           
            alb_c = 0.472+2.*(alphdi-0.472)*(hgbq-1.)       
          else if (hgbq.gt.0.05.and.hgbq.le.1.) then          
                 alb_c = 0.2467+                                        
     &                   (0.7049*hgbq)-(0.8608*(hgbq*hgbq))+
     &                   (0.3812*(hgbq*hgbq*hgbq))                  
          else                                                      
            alb_c = 0.1+3.6*hgbq                                
          endif                                                     
        else                                                        
      ! Melting ice.
      !--------------
          if (hgbq.gt.1.5) then                                  
            alb_c = albice                                           
          else if (hgbq.gt.1..and.hgbq.le.1.5)  then          
                 alb_c = 0.472+(2.*(albice-0.472)*(hgbq-1.))     
          else if (hgbq.gt.0.05.and.hgbq.le.1.) then          
                 alb_c = 0.2467+0.7049*hgbq                          
     &                  -(0.8608*(hgbq*hgbq))
     &                  +(0.3812*(hgbq*hgbq*hgbq)) 
          else                                                      
            alb_c = 0.1+3.6*hgbq
          endif                                                     
        endif                                                       
      endif                                                         
     
      !--------------------------
      ! Correction due to clouds
      !--------------------------
      alb_o= alb_c + cgren                                           

      WRITE(84,*) ' alb_o: ', alb_o
      WRITE(84,*) ' alb_c: ', alb_c

!------------------------------------------------------------------------------!
!- Fin de la routine shine -
      RETURN
      END                                                               
