data(log_om): 'meuse.eas', x=1, y=2, v=4, min=20, max=40, radius=1500, X=4, average;
variogram(log_om): 0.0602857258829195 Nug(0) + 0.207344805460287 Exp(495.103941013073,0,0,0,1,1);
mask: 'dist';
method: gs;
set nsim=50;
predictions(log_om): 'var1.pred.hdr';
