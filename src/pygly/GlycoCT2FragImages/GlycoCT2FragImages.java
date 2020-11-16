
import java.io.File;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.IOException;
import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.lang.RuntimeException;
import java.lang.IllegalArgumentException;

import org.glycoinfo.application.glycanbuilder.converterWURCS2.WURCS2Parser;

import org.eurocarbdb.application.glycanbuilder.Glycan;
import org.eurocarbdb.application.glycanbuilder.Residue;
import org.eurocarbdb.application.glycanbuilder.converterGlycoCT.GlycoCTCondensedParser;
import org.eurocarbdb.application.glycanbuilder.util.GraphicOptions;
import org.eurocarbdb.application.glycanbuilder.massutil.MassOptions;
import org.eurocarbdb.application.glycanbuilder.renderutil.SVGUtils;
import org.eurocarbdb.application.glycanbuilder.renderutil.GlycanRendererAWT;
import org.eurocarbdb.application.glycanbuilder.linkage.Union;
import org.eurocarbdb.application.glycanbuilder.BuilderWorkspace;
import org.eurocarbdb.application.glycanbuilder.Fragmenter;
import org.eurocarbdb.application.glycanbuilder.FragmentCollection;
import org.eurocarbdb.application.glycanbuilder.FragmentEntry;
import org.eurocarbdb.application.glycanbuilder.logutility.LogUtils;


import org.apache.log4j.BasicConfigurator;

public class GlycoCT2FragImages
{
    private static String readFileAsString(String filePath) throws IOException{
	byte[] buffer;
	int maxbuf=10*1024; // 10K, big enough?
	if (filePath.equals("-")) {
	    buffer = new byte[maxbuf];
	} else {
	    buffer = new byte[(int) new File(filePath).length()];
	}
	InputStream f = null;
	int readlen = 0;
	try {
	    if (filePath.equals("-")) {
		f = System.in;
	    } else {
		f = new BufferedInputStream(new FileInputStream(filePath));
	    }
	    readlen = f.read(buffer);
	} finally {
	    if (f != null) try { f.close(); } catch (IOException ignored) { }
	}
	if (readlen >= maxbuf) {
	    throw new IOException("GlycoCTCondensed glycan on standard input too big!");
	}
	return new String(buffer,0,readlen);
    }
    
    private static String changeExtn(String filePath, String extn) throws IOException {
	return changeExtn(filePath,extn,null,true);
    }

    private static String changeExtn(String filePath, String extn, String old) throws IOException {
	return changeExtn(filePath,extn,old,true);
    }

    private static String changeExtn(String filePath, String newextn, 
				     String oldextn, Boolean stripPath) throws IOException {
	if (filePath.equals("-")) {
	    throw new IOException("Can't change extension for standard input");
	}
	File f = new File(filePath);
	String name;
	if (stripPath) {
	    name = f.getName();
	} else {
	    name = f.getPath();
	}
	String extn="";
	String base=name;
	int idx = name.lastIndexOf('.');
	if ((idx >= 1) && (idx <= (name.length()-2))) {
	    base = name.substring(0,idx);
	    extn = name.substring(idx+1);
	}
	if ((oldextn != null) && (extn != oldextn)) {
	    throw new IOException("Bad extension");
	}
	if (newextn == null) {
	    return base;
        }
	return base + "." + newextn;
    }

    private static void setOrientation(BuilderWorkspace t_gwb, String orientation) {
	if (orientation.equals("RL")) {
	    t_gwb.getGraphicOptions().ORIENTATION = GraphicOptions.RL; 
	} else if (orientation.equals("LR")) {
	    t_gwb.getGraphicOptions().ORIENTATION = GraphicOptions.LR; 
	} else if (orientation.equals("BT")) {
	    t_gwb.getGraphicOptions().ORIENTATION = GraphicOptions.BT; 
	} else if (orientation.equals("TB")) {
	    t_gwb.getGraphicOptions().ORIENTATION = GraphicOptions.TB; 
	} else {
	    throw new IllegalArgumentException("Bad orientation option");
	}
    }

    private static void setDisplay(BuilderWorkspace t_gwb, String display) {
	if (display.equals("normalinfo")) {
	    t_gwb.setDisplay(GraphicOptions.DISPLAY_NORMALINFO);
	} else if (display.equals("normal")) {
	    t_gwb.setDisplay(GraphicOptions.DISPLAY_NORMAL);
	} else if (display.equals("compact")) {
	    t_gwb.setDisplay(GraphicOptions.DISPLAY_COMPACT);
	} else {
	    throw new IllegalArgumentException("Bad display option");
	}
    }

    private static void setFragIon(Fragmenter fragger, String ion) {
	if (ion.equals("B")) {
	    fragger.setComputeBFragments(true);
	    fragger.setComputeYFragments(false);
        } else if (ion.equals("b")) {
            fragger.setComputeBFragments(true);
	    fragger.setComputeYFragments(false);
        } else if (ion.equals("Y")) {
	    fragger.setComputeBFragments(false);
            fragger.setComputeYFragments(true);
        } else if (ion.equals("y")) {
	    fragger.setComputeBFragments(false);
            fragger.setComputeYFragments(true);
        } else if (ion.equals("BY")) {
            fragger.setComputeBFragments(true);
            fragger.setComputeYFragments(true);
        } else if (ion.equals("by")) {
            fragger.setComputeBFragments(true);
            fragger.setComputeYFragments(true);
	} else {
	    throw new IllegalArgumentException("Bad ion option");
	}
    }

    private static void setNotation(BuilderWorkspace t_gwb, String notation) {
	if (notation.equals("cfg")) {
	    t_gwb.setNotation(GraphicOptions.NOTATION_CFG);
	} else if (notation.equals("snfg")) {
        t_gwb.setNotation(GraphicOptions.NOTATION_SNFG);
    } else if (notation.equals("cfgbw")) {
	    t_gwb.setNotation(GraphicOptions.NOTATION_CFGBW);
	} else if (notation.equals("cfglink")) {
	    t_gwb.setNotation(GraphicOptions.NOTATION_CFGLINK);
	} else if (notation.equals("uoxf")) {
	    t_gwb.setNotation(GraphicOptions.NOTATION_UOXF);
	} else if (notation.equals("text")) {
	    t_gwb.setNotation(GraphicOptions.NOTATION_TEXT);
	} else if (notation.equals("uoxfcol")) {
	    t_gwb.setNotation(GraphicOptions.NOTATION_UOXFCOL);
	} else {
	    throw new IllegalArgumentException("Bad notation option");
	}
    }

    public static void main(String[] args) throws Exception 
    {
        BasicConfigurator.configure();
	// LogUtils.setGraphicalReport(false);
	BuilderWorkspace t_gwb = new BuilderWorkspace(new GlycanRendererAWT());

	// Single cleavage B/Y ions only
	Fragmenter fragger = new Fragmenter();
	fragger.setComputeAFragments(false);
	fragger.setComputeBFragments(false);
	fragger.setComputeCFragments(false);
	fragger.setComputeXFragments(false);
	fragger.setComputeYFragments(false);
	fragger.setComputeZFragments(false);
	fragger.setComputeInternalFragments(false);
	fragger.setMaxNoCleavages(1);
	fragger.setMaxNoCrossRings(0);
	
	GlycoCTCondensedParser parser = new GlycoCTCondensedParser(true);
	WURCS2Parser wparser = new WURCS2Parser();
	MassOptions mo = new MassOptions();
	boolean mass_opts=false;
	
	setFragIon(fragger,"Y");
	setNotation(t_gwb,"cfg");
	setDisplay(t_gwb,"normalinfo");
	setOrientation(t_gwb,"RL");
	String imagefmt = "png";
	double scale=4.0;
	boolean reducing_end=true;
	boolean opaque=true;
	String outFile = "";

	for (int i=0; i<args.length; i+=1) {

	    if (args[i].equals("format") && args.length > (i+1)) {
		imagefmt = args[i+1];
		i += 1;
		continue;
	    }
	    if (args[i].equals("scale") && args.length > (i+1)) {
		scale = Double.parseDouble(args[i+1]);
		i += 1;
		continue;
	    }
	    if (args[i].equals("redend") && args.length > (i+1)) {
		reducing_end = Boolean.parseBoolean(args[i+1]);
		i += 1;
		continue;
	    }
	    if (args[i].equals("orient") && args.length > (i+1)) {
		setOrientation(t_gwb,args[i+1]);
		i += 1;
		continue;
	    }
	    if (args[i].equals("notation") && args.length > (i+1)) {
		setNotation(t_gwb,args[i+1]);
		i += 1;
		continue;
	    }
	    if (args[i].equals("display") && args.length > (i+1)) {
		setDisplay(t_gwb,args[i+1]);
		i += 1;
		continue;
	    }
	    if (args[i].equals("opaque") && args.length > (i+1)) {
            opaque = Boolean.parseBoolean(args[i+1]);
            i += 1;
            continue;
        }
	    if (args[i].equals("ion") && args.length > (i+1)) {
		setFragIon(fragger,args[i+1]);
		i += 1;
		continue;
	    }
	    if (args[i].equals("base") && args.length > (i+1)) {
		outFile = args[i+1];
		i += 1;
		continue;
	    }

	    String glycanstr = readFileAsString(args[i]);
	    if (outFile.equals("")) {
		outFile = changeExtn(args[i],null);
	    }
	    
	    Glycan glycan;
			    if (glycanstr.startsWith("WURCS")) {
			        glycan = wparser.readGlycan(glycanstr, mo);
			    } else if (glycanstr.startsWith("RES")) {
			        glycan = parser.readGlycan(glycanstr, mo);
			    } else {
			        throw new IllegalArgumentException("Bad glycan descriptor!");
			        //glycan = gparser.readGlycan(glycanstr, mo);
			    }

	    FragmentCollection collection = fragger.computeAllFragments(glycan);
	    int bindex=1;
	    int yindex=1;
	    String name;
	    for (FragmentEntry fe: collection.getFragments()) {
		if (fe.getName().equals("B")) {
		    name = "B" + bindex;
		    bindex ++;
		} else if (fe.getName().equals("Y")) {
		    name = "Y" + yindex;
		    yindex ++;
		} else {
		    continue;
        }

        // SVGUtils.export(t_gwb.getGlycanRenderer(), outFile+'-'+name+'.'+imagefmt, new Union<Glycan>(fe.getFragment()), mass_opts,reducing_end,scale,imagefmt);
	    BufferedImage img = t_gwb.getGlycanRenderer().getImage(glycan, opaque, mass_opts, reducing_end, scale);
        File outputfile = new File(outFile+'-'+name+'.'+imagefmt);
        ImageIO.write(img, imagefmt, outputfile);

		Glycan g = new Glycan(fe.getFragment().getRoot().firstChild(),true,mo);
		String gct = fe.getFragment().toString();
		OutputStream f = null;
		try {
                    f = new BufferedOutputStream(new FileOutputStream(outFile+'-'+name+".txt"));
		    f.write(gct.getBytes());
		    f.close();
	        }
                finally {
		    if (f != null) try { f.close(); } catch (IOException ignored) { }
                }
            }
	    outFile = "";
	}
    }
}
