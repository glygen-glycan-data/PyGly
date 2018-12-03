/* Any JavaScript here will be loaded for all users on every page load. */

// Common action for all pages
window.$ = jQuery;
window.mw = window.mediaWiki = mw;
console.log("Common.js is running");

// s stands for the scripts, while p stand for the loading status for the page
var syncStatus = {"s1": false, "s2": false, "s3": false, "s4": false, "p": false};
function injectEssentialElementForTheViewer(){
    /* Inject the 3 script and the style sheet into each page */

    jQuery.getScript("https://cdn.jsdelivr.net/gh/Hendricks27/Glycan_hierarchical_relationship_viewer/hgv.js", function () {
        syncStatus.s1 = true;
        syncAndContinue();
    });
    jQuery.getScript("https://cdnjs.cloudflare.com/ajax/libs/vis/4.19.1/vis.min.js", function () {
        syncStatus.s2 = true;
        syncAndContinue();
    });
    jQuery.getScript("https://d3js.org/d3.v3.min.js", function () {
        syncStatus.s3 = true;
        syncAndContinue();
    });
    // https://edwardslab.bmcb.georgetown.edu/glycomotifdev/index.php?title=MediaWiki:Components.js&action=raw&ctype=text/javascript
    jQuery.getScript("https://cdn.jsdelivr.net/gh/Hendricks27/presentation/glycomotif/components6.js", function () {
        syncStatus.s4 = true;
        syncAndContinue();
    });
    $('head').append('<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/vis/4.19.1/vis.min.css" type="text/css" />');
}
injectEssentialElementForTheViewer();

jQuery(document).ready(function() {
    syncStatus.p = true;
    syncAndContinue();
});

var hgv_option_template = {
    essentials : {
        div_ID : "", // the ID of div container
        component : {}, // the data, it will be added momentarily
        topoOnly : 0,
        viewRoot : "",
        useGlyTouCanAsImageSource : true,
        GlyTouCanImagePara: {
            style: "extended", // Other Options: normal, compact
            format: "png", // Other Options: jpg
            notation: "cfg" // Other Options: cfgbw, uoxf, uoxf-color, cfg-uoxf, iupac
        },
        imgURL1 : "img/", // Unnecessary if useGlyTouCanAsImageSource is true
        imgURL2 : ".png"
    },
    display : {
        enableTitle: false,
        enableNavi: true,
        naviOption: {
            size: 0.2,
            position: 4
        },
        orientation: 2 // 1, 2, 3, 4 Stand for top2bottom left2right bottom2top right2left
    },
    contextMenu : {
        enable: true,
        defaultMenu: false,
        externalURL1: "", //"https://edwardslab.bmcb.georgetown.edu/glycomotif/GM.",
        externalURL2: ""
    }
};
var divid1 = "viewer_topology";
var divid2 = "viewer_topology_navigator2";
var divid3 = "viewer_topology_navigator";

var pageCollection = document.getElementById("information_for_commonjs").getAttribute("data-collection");
var gtcacc = document.getElementById("information_for_commonjs").getAttribute("data-glytoucan");
var appPrefix = document.getElementById("information_for_commonjs").getAttribute("data-prefix");

hgv_option_template.contextMenu.externalURL1 = "https://edwardslab.bmcb.georgetown.edu/" + appPrefix + "/GM.";

var option1 = jQuery.extend(true, {}, hgv_option_template);
var option2 = jQuery.extend(true, {}, hgv_option_template);
var option3 = jQuery.extend(true, {}, hgv_option_template);
option1.essentials.div_ID = divid1;
option2.essentials.div_ID = divid2;
option3.essentials.div_ID = divid3;
//option1.display.orientation = 2;
option2.essentials.GlyTouCanImagePara.style = "compact";
option3.essentials.GlyTouCanImagePara.style = "compact";

var viewer1, viewer2, viewer3;

function syncAndContinue (){
    var proceed = true;
    for (var k in syncStatus){
        var v = syncStatus[k];
        if (!v){
            proceed = false
        }
    }
    if (proceed){
        console.log("scripts and page loading complete");
        if (pageCollection == "GM"){
            deepCopyViewer();
            fillTheDiv();
        }
    }
}

function deepCopyViewer(){
    viewer1 = jQuery.extend(true, {}, glycanviewer);
    viewer2 = jQuery.extend(true, {}, glycanviewer);
    viewer3 = jQuery.extend(true, {}, glycanviewer);
}

function fillTheDiv() {
    locateViewer1();
    locateViewer2();
    locateViewer3();
}

// Load topology
function locateViewer1(){
    var flag = false;
    for (var key in topologyComponents){
        if (Object.keys(topologyComponents[key].nodes).includes(gtcacc)){
            option1.essentials.component = topologyComponents[key];
            //var style = '<style>#'+ divid1 + '{width: calc(100%); height: calc(70vh); overflow: hidden; border: solid;border-color: lightgrey;}</style>';
            var style = '<style>#'+ divid1 + '{width: calc(100%); height: calc(70vh); overflow: hidden; border: solid;border-color: lightgrey;}</style>';
            $('head').append(style);
            flag = true;
            break;
        }
    }
    if (flag){
        viewer1.init(option1);
        viewer1.network.fit();
        setTimeout(function(){
            viewer1.network.selectNodes([gtcacc]);
        }, 1500);
    }
    else{
        console.log("Bug in load viewer1")
    }
}

// Load sub-structure matching
function locateViewer2(){
    var flag = false;
    for (var key in redonlyComponents){
        if (key == gtcacc){
            option2.essentials.component = redonlyComponents[key];
            option2.essentials.component.nodes[gtcacc].alternativeImageURL = "https://glytoucan.org/glycans/" + gtcacc + "/image?style=extended&format=png&notation=cfg";
            var style = '<style>#'+ divid2 + '{width: calc(100%); height: calc(70vh); overflow: hidden; border: solid;border-color: lightgrey;}</style>';
            $('head').append(style);
            flag = true;
            break;
        }
    }
    if (flag){
        viewer2.init(option2);
        viewer2.network.fit();
        setTimeout(function(){
            viewer2.network.selectNodes([gtcacc]);
        }, 1500);

    }
    else{
        console.log("Bug in load viewer2")
    }
}

function locateViewer3(){
    var flag = false;
    for (var key in nonredonlyComponents){
        if (key == gtcacc){
            option3.essentials.component = nonredonlyComponents[key];
            option3.essentials.component.nodes[gtcacc].alternativeImageURL = "https://glytoucan.org/glycans/" + gtcacc + "/image?style=extended&format=png&notation=cfg";
            var style = '<style>#'+ divid3 + '{width: calc(100%); height: calc(70vh); overflow: hidden; border: solid;border-color: lightgrey;}</style>';
            $('head').append(style);
            flag = true;
            break;
        }
    }
    if (flag){
        viewer3.init(option3);
        viewer3.network.fit();
        setTimeout(function(){
            viewer3.network.selectNodes([gtcacc]);
        }, 1500);

    }
    else{
        console.log("Bug in load viewer3")
    }
}

