# Simulation output sets

| output set    | fAPAR data               |  climate data       | PPFD          |soil moisture   | SM limit.  | temp. ramp  | remarks |
|---------------|--------------------------|---------------------|---------------|----------------|------------|-------------|---------|
| s1            | MOIDS EVI 250 m, monthly | CRU TS 3.22 / WATCH | SPLASH        |SPLASH          | yes        | no          |         |
| s2            | MOIDS EVI 250 m, monthly | CRU TS 3.22 / WATCH | SPLASH        |SPLASH          | no         | no          |         |
| s3            | MOIDS EVI 250 m, monthly | CRU TS 3.22 / WATCH | SPLASH        |SPLASH          | yes        | yes         |         |
| s4            | MOIDS EVI 250 m, monthly | FLUXNET2015 meteo   | SPLASH        |SPLASH          | yes        | yes         |         |
| s5            | MOIDS EVI 250 m, monthly | FLUXNET2015 meteo   | SPLASH        |SPLASH          | no         | yes         |         |
| s6            | MOIDS EVI 250 m, monthly | FLUXNET2015 meteo   | SPLASH        |SWBM            | yes        | yes         | Something went wrong. Check 'Pay' results with internally (SPLASH-) calculated net radiation. |
| s7            | MOIDS EVI 250 m, monthly | FLUXNET2015 meteo   | SPLASH        |SPLASH          | no         | no          |         |
| s8            | MOIDS EVI 250 m, monthly | FLUXNET2015 meteo   | SPLASH        |SPLASH          | yes        | no          |         |
| s9            | MOIDS EVI 250 m, monthly | FLUXNET2015 meteo   | FLUXNET2015   |SPLASH          | no         | no          |         |
| s10           | MOIDS EVI 250 m, daily   | FLUXNET2015 meteo   | FLUXNET2015   |SPLASH          | no         | no          |         |
| s11           | MOIDS EVI 250 m, daily   | FLUXNET2015 meteo   | FLUXNET2015   |SPLASH          | no         | no          | should be identical to s10. re-done climate and fapar input files after FLUXNET2015 update. |
| s12           | MOIDS EVI 250 m, daily   | FLUXNET2015 meteo   | FLUXNET2015   |SWBM            | no         | no          |         |
