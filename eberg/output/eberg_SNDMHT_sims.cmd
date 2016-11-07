data(SNDMHT.t): 'eberg.eas', x=1, y=2, z=3, v=4, min=40, max=50, radius=2500, X=5&6&7&8&9&10&11&12&13&14, average;
variogram(SNDMHT.t): 0 Nug(0) + 0.709678483417632 Exp(128.002672403087,0,0,0,1,4e-04);
set zmap=-0.05;
mask: 'PC1','PC2','PC4','PC5','PC6','PC7','PC8','PC9','PC10','PC11';
predictions(SNDMHT.t): 'var1.pred.hdr';
variances(SNDMHT.t): 'var1.svar.hdr';
