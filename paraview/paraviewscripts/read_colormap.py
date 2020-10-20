import xml.dom.minidom

def colormap(filein):
    
    map=xml.dom.minidom.parse(filein)
    rlist=map.getElementsByTagName('Point')
    rlist2=map.getElementsByTagName('ColorMap')
    colorspace=str(rlist2[0].attributes["space"].value)

    piecewise_fun = []
    lut = []

    for i in range(len(rlist)):

        x=float(rlist[i].attributes["x"].value)
        o=float(rlist[i].attributes["o"].value)
        r=float(rlist[i].attributes["r"].value)
        g=float(rlist[i].attributes["g"].value)
        b=float(rlist[i].attributes["b"].value)
       
        piecewise_fun.append(x)
        piecewise_fun.append(o)

        lut.append(x)
        lut.append(r)       
        lut.append(g)     
        lut.append(b)

    return colorspace, piecewise_fun, lut
