import os

def html_tree(path,title):

    for path,dirs,files in os.walk(path.split("/")[0]):

        if os.path.isdir(path):
            with open(path+"/index.html", 'w') as f:
                f.write("<html>\n <head>\n")
                f.write("  <title>%s</title>\n"%title)
                f.write(" </head>\n")
                f.write(" <body>\n")
                f.write(title)
                f.write("   <br>")
                f.write("   <ul>")
                for d in os.listdir(path):
                    if(os.path.isdir(path+"/"+d)):
                        s = '\t<li><a href="%s/index.html">%s</a></li>\n'%(d,d)
                        f.write(s)
                    elif os.path.isfile(path+"/"+d) and not d.endswith("html"):
                        s = ('\t<li><a href="%s">%s</a></li>'%(d,d))
                        f.write(s)
                f.write("  </ul>")
                f.write("  <br>")
                f.write(" </body>\n")
                f.write("</html>\n")

if __name__ =="__main__":
	html_tree("","List")	
