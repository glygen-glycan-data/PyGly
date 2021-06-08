
import java.io.File;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.FileInputStream;
import java.io.OutputStream;
import java.io.FileOutputStream;
import java.io.BufferedInputStream;
import java.io.IOException;
import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.lang.RuntimeException;
import java.lang.IllegalArgumentException;
import java.lang.StringIndexOutOfBoundsException;
import java.lang.NullPointerException;
import java.lang.StackOverflowError;

import org.glycoinfo.application.glycanbuilder.converterWURCS2.WURCS2Parser;
import org.glycoinfo.GlycanFormatconverter.Glycan.GlycanException;
import org.glycoinfo.WURCSFramework.util.oldUtil.ConverterExchangeException;
import org.glycoinfo.WURCSFramework.util.array.WURCSFormatException;

import org.eurocarbdb.application.glycanbuilder.Glycan;
import org.eurocarbdb.application.glycanbuilder.Residue;
import org.eurocarbdb.application.glycanbuilder.converterGlycoCT.GlycoCTCondensedParser;
import org.eurocarbdb.application.glycanbuilder.util.GraphicOptions;
import org.eurocarbdb.application.glycanbuilder.massutil.MassOptions;
import org.eurocarbdb.application.glycanbuilder.renderutil.SVGUtils;
import org.eurocarbdb.application.glycanbuilder.renderutil.GlycanRendererAWT;
import org.eurocarbdb.application.glycanbuilder.linkage.Union;
import org.eurocarbdb.application.glycanbuilder.BuilderWorkspace;
import org.eurocarbdb.resourcesdb.monosaccharide.MonosaccharideException;
import org.eurocarbdb.MolecularFramework.util.visitor.GlycoVisitorException;


public class GlycoCT2Image
{

	private static String readFileAsString(String filePath) throws IOException{
		byte[] buffer;
		int maxbuf = 10*1024; // 10K, big enough?
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
		// GlycanWorkspace -> BuilderWorkspace: different constructor
		GlycanRendererAWT t_grawt = new GlycanRendererAWT();
		BuilderWorkspace t_gwb = new BuilderWorkspace(t_grawt);

		GlycoCTCondensedParser parser = new GlycoCTCondensedParser(true);
		WURCS2Parser wparser = new WURCS2Parser();

		MassOptions mo = new MassOptions();
		boolean mass_opts=false;

		setNotation(t_gwb,"cfg");
		setDisplay(t_gwb,"normalinfo");
		setOrientation(t_gwb,"RL");
		String imagefmt = "png";
		double scale=4.0;
		boolean reducing_end=true;
		boolean opaque=true;
		boolean force=false;
		String outDir = "";
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
			if (args[i].equals("force") && args.length > (i+1)) {
				force = Boolean.parseBoolean(args[i+1]);
				i += 1;
				continue;
			}
			if (args[i].equals("outdir") && args.length > (i+1)) {
				outDir = args[i+1];
				i += 1;
				continue;
			}
			if (args[i].equals("out") && args.length > (i+1)) {
				outFile = args[i+1];
				i += 1;
				continue;
			}

			String glycanstr = readFileAsString(args[i]);
			if (outFile.equals("")) {
			    if (outDir.equals("")) {
				outFile = changeExtn(args[i],imagefmt);
			    } else {
                                File f = new File(args[i]);
			        String name = f.getName();
			        String newname = changeExtn(name,imagefmt);
				outFile = outDir + File.separator + newname;
			    }
			}

			try {
                            File outputfile = new File(outFile);
			    if (force || !outputfile.exists()) {

			    Glycan glycan;
			    if (glycanstr.startsWith("WURCS")) {
			        glycan = wparser.readGlycan(glycanstr, mo);
			    } else if (glycanstr.startsWith("RES")) {
			        glycan = parser.readGlycan(glycanstr, mo);
			    } else {
			        throw new IllegalArgumentException("Bad glycan descriptor!");
			    }

			    if (imagefmt.equalsIgnoreCase("png") || imagefmt.equalsIgnoreCase("jpg") || imagefmt.equalsIgnoreCase("jpeg")) {
			        BufferedImage img = t_gwb.getGlycanRenderer().getImage(glycan, opaque, mass_opts, reducing_end, scale);
                                ImageIO.write(img, imagefmt, outputfile);
				System.out.println(args[i]);
			    }
			    else if (imagefmt.equalsIgnoreCase("svg")) {

                    String t_svg = SVGUtils.getVectorGraphics(t_grawt, new Union<Glycan>(glycan), mass_opts, reducing_end);

                    FileWriter outputfilewriter = new FileWriter(outFile);
                    outputfilewriter.write(t_svg);
                    outputfilewriter.close();
				System.out.println(args[i]);
			    }
			    else {
			        throw new IllegalArgumentException("Image format " + imagefmt + " is not supported");
			    }
			    }

			}
			catch (StackOverflowError | 
                               WURCSFormatException | 
                               GlycoVisitorException | 
                               MonosaccharideException | 
			       ConverterExchangeException | 
                               GlycanException ex) {
				System.out.println(args[i] + ": " + ex.getClass().getSimpleName() + ": " + ex.getMessage());
			}

			outFile = "";

		}
	}
}
