# trackO2
development of track model for O2

dimanche 20 décembre 2015, 00:51:41 (UTC+0100)
1) Base tracking/updating routines done.
2) test.C is comparison between AliExternalTrackParam and TrackO2.
creates then propagates track to X=400 in steps of 2cm (with Bz only), 
periodically rotating to pt frame, correcting for material and 
updating with measurement, then propagate to x=0 using 3D mag.field
To run it:
aliroot
root [0] gSystem->Load("libTrackO2.so");
root [1] .x test.C+
