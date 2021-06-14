
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
import org.eurocarbdb.resourcesdb.monosaccharide.MonosaccharideException;
import org.eurocarbdb.MolecularFramework.util.visitor.GlycoVisitorException;

public class GlycoCT2GWB
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

	public static void main(String[] args) throws Exception
	{
		// GlycanWorkspace -> BuilderWorkspace: different constructor
		GlycoCTCondensedParser parser = new GlycoCTCondensedParser(true);
		WURCS2Parser wparser = new WURCS2Parser();

		String outFile = "";

		MassOptions mo = new MassOptions();

		for (int i=0; i<args.length; i+=1) {

			String glycanstr = readFileAsString(args[i]);

			try {
			    Glycan glycan;
			    if (glycanstr.startsWith("WURCS")) {
			        glycan = wparser.readGlycan(glycanstr, mo);
			    } else if (glycanstr.startsWith("RES")) {
			        glycan = parser.readGlycan(glycanstr, mo);
			    } else {
			        throw new IllegalArgumentException("Bad glycan descriptor!");
			    }

                    String gwb = glycan.toString();
                    System.out.println(gwb);

			}

			catch (GlycoVisitorException ex) {
				// pass...
			}
			catch (MonosaccharideException ex) {
				// pass...
			}
			catch (RuntimeException ex) {
				// pass...
			}

			outFile = "";

		}
	}
}
