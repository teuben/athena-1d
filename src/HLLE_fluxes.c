#if 0
         roe_average(wl[i],wr[i],b1[i],droe,v1roe,v2roe,v3roe,hroe,b1roe,b2roe,
                     b3roe,x,y,pbl,pbr,el,er);

#if defined(ISOTHERMAL) AND defined(HYDRO)
         eigen_iso_hydro(v1roe,v2roe,v3roe,ev,rem,lem);
#endif
#if defined(ADIABATIC) AND defined(HYDRO)
         adiabatic_hydro_conserved(u_roe,ev,rem,lem);
#endif
#if defined(ISOTHERMAL) AND defined(MHD)
         eigen_iso_mhd(u_roe,b1,x,y,'C',ev,rem,lem);
#endif
#if defined(ADIABATIC) AND defined(MHD)
         eigen_ad_mhd(u_roe,b1,x,y,'C',ev,rem,lem);
#endif
         maxevroe = MAX(maxevroe,(MAX(fabs(ev[0]),fabs(ev[NVAR-1]))));
      }

/* Construct left- and right-fluxes from left- and right-state vectors */
/* (eq. XX)  */

      lr_fluxes();


      cl = sqrt((GAMMA*pl + 2.0*pbl)/dl);
      cr = sqrt((GAMMA*pr + 2.0*pbr)/dr);
      bp = MAX(MAX(ev[NVAR-1],(v1r + cr)), 0.0);
      bm = MIN(MIN(ev[0]     ,(v1l - cl)), 0.0);
      for (n=0; n<NVAR; n++) {
         f[i][n] = ((bp*fl[n]-bm*fr[n]) + bp*bm*(ur[i][n]-ul[i][n]))/(bp-bm);
      }
   }
   return(maxevroe);
}
#endif
