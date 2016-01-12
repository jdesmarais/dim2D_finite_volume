import xml.etree.ElementTree as ET

if __name__=="__main__":
    
    tree = ET.parse('file.xml')
    root = tree.getroot()

    for window in root.findall('window'):
        
        for plot in window.findall('plot'):
            print plot.find('var').text
    
