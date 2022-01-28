from ROOT import *
from array import array
outfile = TFile("Output_electron_trigger_sf.root","RECREATE")
i=1
xbins = array('d',[-2.5,-1.566,-1.4442,0,1.4442,1.556,2.5])
ybins = array('d',[0,50,55,60,65,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,300,500,2000])
hist_18 = TH2D("2018","2018",len(xbins)-1,xbins,len(ybins)-1,ybins)

hist_18.SetBinContent(1,i+1,0.981798)
hist_18.SetBinContent(3,i+1,0.961154)
hist_18.SetBinContent(4,i+1,0.953125)
hist_18.SetBinContent(6,i+1,0.980549)

hist_18.SetBinContent(1,i+2,0.977528)
hist_18.SetBinContent(3,i+2,0.964796)
hist_18.SetBinContent(4,i+2,0.954646)
hist_18.SetBinContent(6,i+2,0.976325)

hist_18.SetBinContent(1,i+3,0.981069)
hist_18.SetBinContent(3,i+3,0.96179)
hist_18.SetBinContent(4,i+3,0.956044)
hist_18.SetBinContent(6,i+3,0.986486)

hist_18.SetBinContent(1,i+4,0.98234)
hist_18.SetBinContent(3,i+4,0.966487)
hist_18.SetBinContent(4,i+4,0.954048)
hist_18.SetBinContent(6,i+4,0.995531)

hist_18.SetBinContent(1,i+5,0.984683)
hist_18.SetBinContent(3,i+5,0.970842)
hist_18.SetBinContent(4,i+5,0.961957)
hist_18.SetBinContent(6,i+5,0.996681)

hist_18.SetBinContent(1,i+6,0.981622)
hist_18.SetBinContent(3,i+6,0.973176)
hist_18.SetBinContent(4,i+6,0.96129)
hist_18.SetBinContent(6,i+6,0.994493)

hist_18.SetBinContent(1,i+7,0.984848)
hist_18.SetBinContent(3,i+7,0.976521)
hist_18.SetBinContent(4,i+7,0.969893)
hist_18.SetBinContent(6,i+7,0.988069)

hist_18.SetBinContent(1,i+8,0.975322)
hist_18.SetBinContent(3,i+8,0.969182)
hist_18.SetBinContent(4,i+8,0.962766)
hist_18.SetBinContent(6,i+8,0.993457)

hist_18.SetBinContent(1,i+9,0.968912)
hist_18.SetBinContent(3,i+9,0.972746)
hist_18.SetBinContent(4,i+9,0.959077)
hist_18.SetBinContent(6,i+9,0.980707)

hist_18.SetBinContent(1,i+10,0.978831)
hist_18.SetBinContent(3,i+10,0.98471)
hist_18.SetBinContent(4,i+10,0.977597)
hist_18.SetBinContent(6,i+10,0.97762)

hist_18.SetBinContent(1,i+11,0.995964)
hist_18.SetBinContent(3,i+11,0.992886)
hist_18.SetBinContent(4,i+11,0.98073)
hist_18.SetBinContent(6,i+11,0.992915)

hist_18.SetBinContent(1,i+12,0.99196)
hist_18.SetBinContent(3,i+12,0.992901)
hist_18.SetBinContent(4,i+12,0.976744)
hist_18.SetBinContent(6,i+12,0.992944)

hist_18.SetBinContent(1,i+13,0.991944)
hist_18.SetBinContent(3,i+13,0.983855)
hist_18.SetBinContent(4,i+13,0.988832)
hist_18.SetBinContent(6,i+13,0.981)

hist_18.SetBinContent(1,i+14,0.996)
hist_18.SetBinContent(3,i+14,0.985844)
hist_18.SetBinContent(4,i+14,0.980827)
hist_18.SetBinContent(6,i+14,0.983)

hist_18.SetBinContent(1,i+15,0.99798)
hist_18.SetBinContent(3,i+15,0.991886)
hist_18.SetBinContent(4,i+15,0.98583)
hist_18.SetBinContent(6,i+15,0.989)

hist_18.SetBinContent(1,i+16,0.985)
hist_18.SetBinContent(3,i+16,0.996973)
hist_18.SetBinContent(4,i+16,0.971803)
hist_18.SetBinContent(6,i+16,0.978)

hist_18.SetBinContent(1,i+17,0.997984)
hist_18.SetBinContent(3,i+17,0.998981)
hist_18.SetBinContent(4,i+17,0.976861)
hist_18.SetBinContent(6,i+17,1.00302)

hist_18.SetBinContent(1,i+18,0.995)
hist_18.SetBinContent(3,i+18,0.990991)
hist_18.SetBinContent(4,i+18,0.988956)
hist_18.SetBinContent(6,i+18,0.991)

hist_18.SetBinContent(1,i+19,0.998)
hist_18.SetBinContent(3,i+19,0.998)
hist_18.SetBinContent(4,i+19,1.001)
hist_18.SetBinContent(6,i+19,0.992)

hist_18.SetBinContent(1,i+20,0.996)
hist_18.SetBinContent(3,i+20,0.999)
hist_18.SetBinContent(4,i+20,0.997)
hist_18.SetBinContent(6,i+20,0.998)

hist_18.SetBinContent(1,i+21,0.997)
hist_18.SetBinContent(3,i+21,0.999)
hist_18.SetBinContent(4,i+21,1.00403)
hist_18.SetBinContent(6,i+21,0.994)

hist_18.SetBinContent(1,i+22,0.997)
hist_18.SetBinContent(3,i+22,0.997)
hist_18.SetBinContent(4,i+22,0.998)
hist_18.SetBinContent(6,i+22,0.997)

hist_18.SetBinContent(1,i+23,0.999)
hist_18.SetBinContent(3,i+23,0.999)
hist_18.SetBinContent(4,i+23,0.993)
hist_18.SetBinContent(6,i+23,0.999)


hist_18.SetBinContent(1,i+24,0.998)
hist_18.SetBinContent(3,i+24,0.997)
hist_18.SetBinContent(4,i+24,0.999)
hist_18.SetBinContent(6,i+24,0.998)

hist_18.SetBinContent(1,i+25,0.998)
hist_18.SetBinContent(3,i+25,0.996)
hist_18.SetBinContent(4,i+25,0.996)
hist_18.SetBinContent(6,i+25,0.998)
hist_18.Write()
hist_17 = TH2D("2017","2017",len(xbins)-1,xbins,len(ybins)-1,ybins)

hist_17.SetBinContent(1,i+1,0.909195)
hist_17.SetBinContent(3,i+1,0.948315)
hist_17.SetBinContent(4,i+1,0.937642)
hist_17.SetBinContent(6,i+1,0.896287)

hist_17.SetBinContent(1,i+2,0.912162)
hist_17.SetBinContent(3,i+2,0.951002)
hist_17.SetBinContent(4,i+2,0.939258)
hist_17.SetBinContent(6,i+2,0.898534)

hist_17.SetBinContent(1,i+3,0.920314)
hist_17.SetBinContent(3,i+3,0.961154)
hist_17.SetBinContent(4,i+3,0.938753)
hist_17.SetBinContent(6,i+3,0.912752)

hist_17.SetBinContent(1,i+4,0.926009)
hist_17.SetBinContent(3,i+4,0.956092)
hist_17.SetBinContent(4,i+4,0.941372)
hist_17.SetBinContent(6,i+4,0.9129)

hist_17.SetBinContent(1,i+5,0.937847)
hist_17.SetBinContent(3,i+5,0.959474)
hist_17.SetBinContent(4,i+5,0.944201)
hist_17.SetBinContent(6,i+5,0.935912)

hist_17.SetBinContent(1,i+6,0.951542)
hist_17.SetBinContent(3,i+6,0.957082)
hist_17.SetBinContent(4,i+6,0.947939)
hist_17.SetBinContent(6,i+6,0.941886)

hist_17.SetBinContent(1,i+7,0.932682)
hist_17.SetBinContent(3,i+7,0.959402)
hist_17.SetBinContent(4,i+7,0.949134)
hist_17.SetBinContent(6,i+7,0.970787)

hist_17.SetBinContent(1,i+8,0.952747)
hist_17.SetBinContent(3,i+8,0.954207)
hist_17.SetBinContent(4,i+8,0.952637)
hist_17.SetBinContent(6,i+8,0.926486)

hist_17.SetBinContent(1,i+9,0.930159)
hist_17.SetBinContent(3,i+9,0.950682)
hist_17.SetBinContent(4,i+9,0.944152)
hist_17.SetBinContent(6,i+9,0.93038)

hist_17.SetBinContent(1,i+10,0.951807)
hist_17.SetBinContent(3,i+10,0.963489)
hist_17.SetBinContent(4,i+10,0.973388)
hist_17.SetBinContent(6,i+10,0.952573)

hist_17.SetBinContent(1,i+11,0.967)
hist_17.SetBinContent(3,i+11,0.984694)
hist_17.SetBinContent(4,i+11,0.964682)
hist_17.SetBinContent(6,i+11,0.953722)

hist_17.SetBinContent(1,i+12,0.970766)
hist_17.SetBinContent(3,i+12,0.993884)
hist_17.SetBinContent(4,i+12,0.968718)
hist_17.SetBinContent(6,i+12,0.971944)

hist_17.SetBinContent(1,i+13,0.980847)
hist_17.SetBinContent(3,i+13,0.977845)
hist_17.SetBinContent(4,i+13,0.962739)
hist_17.SetBinContent(6,i+13,0.950655)

hist_17.SetBinContent(1,i+14,0.967742)
hist_17.SetBinContent(3,i+14,0.983773)
hist_17.SetBinContent(4,i+14,0.977845)
hist_17.SetBinContent(6,i+14,0.955912)

hist_17.SetBinContent(1,i+15,0.951)
hist_17.SetBinContent(3,i+15,0.984894)
hist_17.SetBinContent(4,i+15,0.973764)
hist_17.SetBinContent(6,i+15,0.962702)

hist_17.SetBinContent(1,i+16,0.983952)
hist_17.SetBinContent(3,i+16,0.974722)
hist_17.SetBinContent(4,i+16,1.00831)
hist_17.SetBinContent(6,i+16,0.976)

hist_17.SetBinContent(1,i+17,0.978809)
hist_17.SetBinContent(3,i+17,0.973)
hist_17.SetBinContent(4,i+17,0.982794)
hist_17.SetBinContent(6,i+17,0.976)

hist_17.SetBinContent(1,i+18,0.992)
hist_17.SetBinContent(3,i+18,0.989)
hist_17.SetBinContent(4,i+18,0.978)
hist_17.SetBinContent(6,i+18,0.988)

hist_17.SetBinContent(1,i+19,0.994)
hist_17.SetBinContent(3,i+19,0.998)
hist_17.SetBinContent(4,i+19,0.997)
hist_17.SetBinContent(6,i+19,0.986)

hist_17.SetBinContent(1,i+20,0.995)
hist_17.SetBinContent(3,i+20,0.999)
hist_17.SetBinContent(4,i+20,0.993)
hist_17.SetBinContent(6,i+20,0.994)

hist_17.SetBinContent(1,i+21,0.994)
hist_17.SetBinContent(3,i+21,0.998)
hist_17.SetBinContent(4,i+21,0.998)
hist_17.SetBinContent(6,i+21,0.995)

hist_17.SetBinContent(1,i+22,0.99)
hist_17.SetBinContent(3,i+22,0.999)
hist_17.SetBinContent(4,i+22,0.995)
hist_17.SetBinContent(6,i+22,0.994)

hist_17.SetBinContent(1,i+23,0.994)
hist_17.SetBinContent(3,i+23,1.)
hist_17.SetBinContent(4,i+23,0.997)
hist_17.SetBinContent(6,i+23,0.998)

hist_17.SetBinContent(1,i+24,0.997)
hist_17.SetBinContent(3,i+24,0.998)
hist_17.SetBinContent(4,i+24,0.998)
hist_17.SetBinContent(6,i+24,0.997)

hist_17.SetBinContent(1,i+25,0.997)
hist_17.SetBinContent(3,i+25,0.993)
hist_17.SetBinContent(4,i+25,0.993)
hist_17.SetBinContent(6,i+25,0.997)

hist_17.Write()
hist_16 = TH2D("2016","2016",len(xbins)-1,xbins,len(ybins)-1,ybins)

hist_16.SetBinContent(1,i+1,0.94133)
hist_16.SetBinContent(3,i+1,0.995261)
hist_16.SetBinContent(4,i+1,0.991726)
hist_16.SetBinContent(6,i+1,0.944301)

hist_16.SetBinContent(1,i+2,0.95)
hist_16.SetBinContent(3,i+2,0.991879)
hist_16.SetBinContent(4,i+2,0.988426)
hist_16.SetBinContent(6,i+2,0.94088)

hist_16.SetBinContent(1,i+3,0.937811)
hist_16.SetBinContent(3,i+3,0.994279)
hist_16.SetBinContent(4,i+3,0.985194)
hist_16.SetBinContent(6,i+3,0.936352)

hist_16.SetBinContent(1,i+4,0.941032)
hist_16.SetBinContent(3,i+4,0.99209)
hist_16.SetBinContent(4,i+4,0.993213)
hist_16.SetBinContent(6,i+4,0.920863)

hist_16.SetBinContent(1,i+5,0.926014)
hist_16.SetBinContent(3,i+5,0.992239)
hist_16.SetBinContent(4,i+5,0.981195)
hist_16.SetBinContent(6,i+5,0.933649)

hist_16.SetBinContent(1,i+6,0.92907)
hist_16.SetBinContent(3,i+6,0.991238)
hist_16.SetBinContent(4,i+6,0.980435)
hist_16.SetBinContent(6,i+6,0.922811)

hist_16.SetBinContent(1,i+7,0.934633)
hist_16.SetBinContent(3,i+7,0.99026)
hist_16.SetBinContent(4,i+7,0.974138)
hist_16.SetBinContent(6,i+7,0.909808)

hist_16.SetBinContent(1,i+8,0.932203)
hist_16.SetBinContent(3,i+8,0.989236)
hist_16.SetBinContent(4,i+8,0.982814)
hist_16.SetBinContent(6,i+8,0.917503)

hist_16.SetBinContent(1,i+9,0.910148)
hist_16.SetBinContent(3,i+9,0.987435)
hist_16.SetBinContent(4,i+9,0.966976)
hist_16.SetBinContent(6,i+9,0.938624)

hist_16.SetBinContent(1,i+10,0.974849)
hist_16.SetBinContent(3,i+10,0.98288)
hist_16.SetBinContent(4,i+10,0.98191)
hist_16.SetBinContent(6,i+10,0.9749)

hist_16.SetBinContent(1,i+11,1.00708)
hist_16.SetBinContent(3,i+11,0.984909)
hist_16.SetBinContent(4,i+11,0.980847)
hist_16.SetBinContent(6,i+11,0.995976)

hist_16.SetBinContent(1,i+12,1.00302)
hist_16.SetBinContent(3,i+12,0.996964)
hist_16.SetBinContent(4,i+12,0.978979)
hist_16.SetBinContent(6,i+12,0.993)

hist_16.SetBinContent(1,i+13,1.01939)
hist_16.SetBinContent(3,i+13,0.986922)
hist_16.SetBinContent(4,i+13,0.980866)
hist_16.SetBinContent(6,i+13,0.991)

hist_16.SetBinContent(1,i+14,0.999)
hist_16.SetBinContent(3,i+14,0.983936)
hist_16.SetBinContent(4,i+14,0.978937)
hist_16.SetBinContent(6,i+14,0.998)

hist_16.SetBinContent(1,i+15,0.998)
hist_16.SetBinContent(3,i+15,0.993958)
hist_16.SetBinContent(4,i+15,0.981964)
hist_16.SetBinContent(6,i+15,0.998)

hist_16.SetBinContent(1,i+16,0.997)
hist_16.SetBinContent(3,i+16,0.993988)
hist_16.SetBinContent(4,i+16,0.997992)
hist_16.SetBinContent(6,i+16,0.997)

hist_16.SetBinContent(1,i+17,0.997)
hist_16.SetBinContent(3,i+17,0.997)
hist_16.SetBinContent(4,i+17,0.997)
hist_16.SetBinContent(6,i+17,0.996)

hist_16.SetBinContent(1,i+18,0.996)
hist_16.SetBinContent(3,i+18,0.999)
hist_16.SetBinContent(4,i+18,1.)
hist_16.SetBinContent(6,i+18,0.996)

hist_16.SetBinContent(1,i+19,0.994)
hist_16.SetBinContent(3,i+19,0.999)
hist_16.SetBinContent(4,i+19,0.99)
hist_16.SetBinContent(6,i+19,0.994)

hist_16.SetBinContent(1,i+20,0.992)
hist_16.SetBinContent(3,i+20,0.996)
hist_16.SetBinContent(4,i+20,0.989)
hist_16.SetBinContent(6,i+20,0.994)

hist_16.SetBinContent(1,i+21,0.991)
hist_16.SetBinContent(3,i+21,0.997)
hist_16.SetBinContent(4,i+21,0.998)
hist_16.SetBinContent(6,i+21,0.987)

hist_16.SetBinContent(1,i+22,0.99)
hist_16.SetBinContent(3,i+22,0.998)
hist_16.SetBinContent(4,i+22,0.993)
hist_16.SetBinContent(6,i+22,0.989)

hist_16.SetBinContent(1,i+23,0.996)
hist_16.SetBinContent(3,i+23,0.994)
hist_16.SetBinContent(4,i+23,0.997)
hist_16.SetBinContent(6,i+23,0.996)

hist_16.SetBinContent(1,i+24,0.995)
hist_16.SetBinContent(3,i+24,0.999)
hist_16.SetBinContent(4,i+24,0.996)
hist_16.SetBinContent(6,i+24,0.995)

hist_16.SetBinContent(1,i+25,0.995)
hist_16.SetBinContent(3,i+25,0.993)
hist_16.SetBinContent(4,i+25,0.993)
hist_16.SetBinContent(6,i+25,0.995)

hist_16.Write()
outfile.Close()
