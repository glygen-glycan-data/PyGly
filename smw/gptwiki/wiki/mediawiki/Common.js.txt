var curl = window.location.href.split("/");
var pagename = curl[curl.length-1];
var paraFromPage = [];

if ( pagename.startsWith("TG") ){
    checkElement();
}

function checkElement() {
    var pot = document.getElementsByClassName("mw-parser-output");
    var msl = pot[0].getElementsByTagName("p");
    var msele = msl[msl.length - 1];

    var lastprecscan = "";
    for (var i in msele.innerText.split("\n")){

        var sth = msele.innerText.split("\n")[i];
        if (sth == ""){
            continue
        }
        var scandetails = sth.split(";");
        var precscan = scandetails[3];
        if (precscan != lastprecscan) {
            paraFromPage.push({title: "Scan: "+precscan, scan: precscan, annotations: false, zoomHeight: true});
            lastprecscan = precscan;
        }
        var title = "";
        if (scandetails.length == 2) {
            title = "Scan: "+scandetails[0]+" ("+scandetails[1]+")";
        } else {
            if (scandetails[3] == "") {
                title = scandetails[1]+" @ "+scandetails[2]+" min, m/z "+scandetails[4];
            } else {
                title = scandetails[1]+" @ "+scandetails[2]+" min, m/z "+scandetails[4]+" (Scan: " + scandetails[0] + ", Precursor scan: " +scandetails[3]+")";
            }
        }
        paraFromPage.push({title: title, scan: scandetails[0], annotations: true});
    }

    msele.innerHTML = "";
}

if (!( Object.keys(paraFromPage).length === 0 && paraFromPage.constructor === Object)){
    injectResources();
}


function injectResources() {

    var external_resources = [
        "https://cdnjs.cloudflare.com/ajax/libs/d3/5.12.0/d3.min.js",
        "https://cdn.jsdelivr.net/gh/glygen-glycan-data/JSWidgets/MS_Viewer/spectrum-parsers.js",
        "https://cdn.jsdelivr.net/gh/glygen-glycan-data/JSWidgets/MS_Viewer/MSV.js",
        "https://cdn.jsdelivr.net/gh/glygen-glycan-data/JSWidgets/MS_Viewer/util.js"
    ];

    function getScript(url) {
        return new Promise(function (resolve, reject){
            jQuery.ajax({
                type: "GET",
                url: url,
                success: function (d) {
                    resolve(d);
                },
                dataType: "script",
                cache: true
            });
        });
    }
    $("head").append("<link rel='stylesheet' href='https://cdn.jsdelivr.net/gh/glygen-glycan-data/JSWidgets/MS_Viewer/spectrum-viewer.css' type='text/css'>");

    var js_pool = [];
    for (var i in external_resources){
        var jsurl = external_resources[i];
        js_pool.push(getScript(jsurl));
    }

    Promise.all(js_pool).then(getScript("https://cdnjs.cloudflare.com/ajax/libs/d3-tip/0.9.1/d3-tip.min.js"))
        .then(function () {
            loadSVG();
        })
}

function loadSVG() {
    // console.log("Start drawing");

    var pot = document.getElementsByClassName("mw-parser-output")[0];

    var msContainer = document.createElement("div");
    msContainer.setAttribute("id", "specpanel");
    pot.appendChild(msContainer);

    var peptideAcc = document.getElementById("msv_para").getAttribute("data-peptide");
    var charger = document.getElementById("msv_para").getAttribute("data-z1");
    var pmz = document.getElementById("msv_para").getAttribute("data-mz1");
    var spectra_folder = document.getElementById("msv_para").getAttribute("data-spectra");
    var showxic = document.getElementById("msv_para").getAttribute("data-showxic") || "true";

    // convert to boolean...
    showxic = (showxic.trim().toLowerCase() != "false");

    var showcycle = document.getElementById("msv_para").getAttribute("data-showcycle") || "false,true,true";

    // convert to array of boolean...
    showcycle = showcycle.split(",");
    var showcyclearray = [];
    for (var i in showcycle) {
        showcyclearray.push((showcycle[i].trim().toLowerCase() != "false"));
    }
    var msmsv = msv();
    var param0 = {
        spectra: "",
        format: "json",
        graphtype: "chromatogram",
        zoomHeight: true,
        width: pot.clientWidth * 0.995,
        titleTag: "h4"
    };

    if (showxic) {
        param0["title"] = ["XIC: m/z "+pmz];
        param0["spectra"] = "https://edwardslab.bmcb.georgetown.edu/~nedwards/dropbox/pBYmLSkGeq/" + spectra_folder + "/" + peptideAcc + "." + charger + '.json';
        msmsv.addLabelledSpectrum('specpanel','chrom', param0);
    }

    for (var i in paraFromPage) {
        var sc = paraFromPage[i]["scan"];
        var title = paraFromPage[i]["title"];
        var params = {
            scan: sc,
            format: "json",
            spectra: "https://edwardslab.bmcb.georgetown.edu/~nedwards/dropbox/pBYmLSkGeq/" + spectra_folder + "/" + sc + '.json',
            width: pot.clientWidth * 0.995,
            title: title,
            titleTag: "h4",
            zoomgroup: 3
        };

        if (paraFromPage[i]["annotations"]) {
            params["annotations"] = "https://edwardslab.bmcb.georgetown.edu/~nedwards/dropbox/pBYmLSkGeq/annotations/" + peptideAcc + "." + charger + ".json";
        }
        if (paraFromPage[i]["zoomHeight"]) {
            params["zoomHeight"] = true;
        } else {
            params["zoomHeight"] = false;
        }

        params["show"] = showcyclearray[i] == true;

        msmsv.addLabelledSpectrum('specpanel','spec'+sc, params);
    }

    msmsv.done();
}