#!/usr/bin/env python

"""
The script reads an sqlite3 database produced by CADEE
and creates an interactive web UI application for
visualizing and analysing the results (index.html).

usage: python analyse.py cadee.db
       this will create your index.html

Author: {0} ({1})

This program is part of CADEE, the framework for
Computer-Aided Directed Evolution of Enzymes.

"""


from __future__ import print_function
import sys
import os
import sqlite3
import json

__author__ = "Miha Purg, Beat Amrein"
__email__ = "miha.purg@gmail.com, beat.amrein@gmail.com"

html="""
<!doctype html>

<html lang="en">
<head>
  <meta charset="utf-8">

  <title>No title</title>
  <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/jqueryui/1.12.0/jquery-ui.css">  
  <link rel="stylesheet" href="http://yui.yahooapis.com/pure/0.6.0/pure-min.css">

  <script type="text/javascript" src="https://cdn.plot.ly/plotly-latest.min.js"></script>
  <script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.1.0/jquery.js"></script>
  <script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/jqueryui/1.12.0/jquery-ui.js"></script>

  <!--[if lt IE 9]>
    <script src="http://html5shiv.googlecode.com/svn/trunk/html5.js"></script>
  <![endif]-->
  <style type="text/css">
* {
    font-family: "arial";
    -webkit-box-sizing: border-box;
       -moz-box-sizing: border-box;
            box-sizing: border-box;
            border-width: 0px;
}

.buttons { 
   margin-left: 20px;
   margin-top: 20px;
}

.sort-button {
   border-radius: 5px;
   background: #88aacc; 
   color: #fafaff;
}

#info-div {
    opacity: 0.9;
    background-color: #fafafa;
    display: none;
    max-width: 500px;
    margin: 0 auto;
    padding: 0.5%;
    border-radius: 5px;
    border: 1px solid #c0c0c0;
    font-size: 12px;
    z-index: 1000;
}

#overlay-div {
    opacity: 0.7;
    background-color: #000;
    width: 100%;
    height: 100%;
    z-index: 500;
    position: absolute;
    left: 0;
    top: 0;
    display: none;

}

#info-header {
    padding: 10px;
    background-color: #fafafa;
    margin-bottom: 10px;
    border-radius: 5px;
    border: 1px solid #d0d0d0;
}

#info-stats {
    padding-left: 10px;
    width: 59%;
    display: inline-block;
    vertical-align: top;
}

#info-actions {
    padding: 0 0 20px 20px;
    width: 40%;
    margin-left: 1%;
    display: inline-block;
    vertical-align: top;
    border-radius: 5px;
    border: 1px solid #d0d0d0;

    background-color: #f6f6f6;
}

input[type=checkbox]
{
  /* 1.5z-sized Checkboxes */
  -ms-transform: scale(1.3); /* IE */
  -moz-transform: scale(1.3); /* FF */
  -webkit-transform: scale(1.3); /* Safari and Chrome */
  -o-transform: scale(1.3); /* Opera */
  transform: scale(1.3);
  margin: 0 10px 10px 0;
}

.actions  {
  text-align: left;
}

input[type="text"] {
  font-family: sans-serif;
  font-size: 12px;
  width: 90%;
  padding: 5px;
}

#output-div {
  margin: 20px;
  padding: 10px;
  border-radius: 5px;
  background-color: #f0f0f0;
}

#output-div > p {
    margin: 10px 0 0 0;
    padding: 10px;
    background-color: #fafafa;
    border-radius: 5px;
    border: 1px solid #c0c0c0;
}

#top_info {
  padding: 10px;
  background-color: #ffaa88;
  opacity: 0.5;
}
</style>


</head>

<body>
<div id="top_info">Info: Use the 'shift' key to toggle the info box, mouse click to select.</div>
<div class="buttons"><button id="sort-barrier" class="sort-button pure-button">Sort by dG#</button>
<button id="sort-name" class="sort-button pure-button">Sort by name</button>
<button id="sort-action" class="sort-button pure-button">Sort by actions</button></div>
<div id="overlay-div"></div>
<div style="overflow: scroll; width:100%;"><div id="graph"></div> </div>
<div id="output-div">CADEE command<p>Nothing yet...</p></div>
<div id="info-div">
    <div id="info-header"></div><div id="info-stats"></div><div id="info-actions"> 
    <p><b>Actions</b></p>
    <form class="actions pure-form" action="javascript:void(0);">
         <label for="inp_SATURATE"><input type="checkbox" name="saturate" id="inp_SATURATE" value="SATURATE"/>Saturate (20AA)</label> <br>          
         <label for="inp_APOLAR"><input type="checkbox" name="apolar" id="inp_APOLAR" value="APOLAR"/>Apolar (8AA)</label>  
         <br>        
         <label for="inp_POLAR"><input type="checkbox" name="polar" id="inp_POLAR" value="POLAR"/>Polar (9AA)</label>
         <br>
         <label for="inp_CUSTOM"><input type="checkbox" name="custom" id="inp_CUSTOM" value="CUSTOM"/>Custom</label> 
         <input type="text" name="CUSTOM_text" id="CUSTOM_text" disabled/>         
    </form>
    </div>
</div>
<div id="comptime-div"></div>
"""

vis="""<script>
var plotted = false;
draw_plot();
var curX;
var curY;
var x_hovered;   // index of element in focus (hovered over - plotly_hover)
var info_visible = false;

data.forEach(function(dp) {
// used to store actions 
    dp.cadee_actions = {
        librarymut: [],
        librarymut_custom: ""
    }
});



function draw_plot() {
    var xs1 = [],
        xs2 = [],
        xs3 = [],
        ys1 = [],
        ys2 = [],
        ys3 = [],
        names_list = []

    data.forEach( function(dp) {
        for (i in dp.barrier) { ys1.push(dp.barrier[i]); xs1.push(dp.name) };
        for (i in dp.exotherm) { ys2.push(dp.exotherm[i]); xs2.push(dp.name) };
        for (i in dp.revbarr) { ys3.push(dp.exotherm[i]); xs3.push(dp.name) };
        names_list.push(dp.name);
    });

    // data1,data2 and data3 are traces for dg#, dg0 and dg#_reverse, respectively
    // each of them contain  number_of_points_per_mutant*number_of_mutant_names  points in 'y' and the same amount of names in 'x' 
    // (x: [name1, name1, name1, name2, name2, name2 ] if 2 names and 3 points per name )
    // (y: [pt1_1, pt1_2, pt1_3, pt2_1, pt2_2, pt2_3]
    // this kind of "group" plotting is enabled with 'boxmode: group' in the layout
    // see: https://plot.ly/javascript/box-plots/
    var data1 = {
            name: 'dG#',
            y: ys1,
            x: xs1,
            type: 'box',
            hoverinfo: "none",
            marker: { color: 'gray' },

            // used for sorting (see sort-barrier function)
            names_list: names_list
    };
    var data2 = {
            name: 'dG0',
            y: ys2,
            x: xs2,
            type: 'box',
            hoverinfo: "none",
            marker: { color: '#ccddcc' }
    };
    var data3 = {
            name: 'dG#r',
            y: ys3,
            x: xs3,
            type: 'box',
            hoverinfo: "none",
            marker: { color: '#ddddff' }
    };


    if (!plotted) {
        //Format the layout
        var layout = {
            xaxis: {
                showgrid: false,
                tickangle: 60,
                showticklabels: true,
            },
            yaxis: {
                title: "Free energy (rel to WT) [kcal/mol]",
                zeroline: true,
                gridcolor: 'rgb(220,220,220)'
            },
              autosize: false,
              width: 40*data.length,
            margin: {
                t: 30,
                b: 150,

            },
            paper_bgcolor: 'rgb(255,255,255)',
            plot_bgcolor: 'rgb(253,253,253)',
            boxmode: "group",
        };

        Plotly.newPlot('graph', [data1, data2, data3], layout);
        plotted = true;
    } else {
        var graph = document.getElementById("graph");
        graph.data = [data1, data2, data3];
        Plotly.redraw(graph);
    };
};

function update_cadee_cmd() {
    var lmuts = [];
    data.forEach(function(dp) {
        if (dp.name == "wt") { return };
        var resid = dp.name.slice(1,-1);

        dp.cadee_actions.librarymut.forEach(function(lm) {
            lmuts.push(resid + ":" + lm); 
        });

        if (dp.cadee_actions.librarymut_custom != "") {
            lmuts.push(resid + ":'" + dp.cadee_actions.librarymut_custom + "'");
        };
    });
    var cmd_str = "Nothing yet..."
    if (lmuts.length > 0) {
        cmd_str = "cadee.py --librarymut " + lmuts.join(" ")
    }
    $("#output-div > p").text(cmd_str);
};

function dock_info_div() {
    var idiv = $("#info-div"); 

    $("#overlay-div").fadeIn("fast");
    idiv.fadeIn("fast");
    idiv.css("opacity", "1")
};

function undock_info_div() {
    $("#overlay-div").hide();
    $("#info-div").css("opacity", "0.9");
};


$("#sort-barrier").click( function() {
    var calcdata = document.getElementById("graph").calcdata[0]; // just the first set - dG#
    //console.log(calcdata);

    var cd_indexes = calcdata.map(function(cd,i) { return i; });

    cd_indexes.sort( function(a,b) { return calcdata[a].med - calcdata[b].med } );
    var graph_data1 = document.getElementById('graph').data[0] // just the first set - dG#
    var cd_names = cd_indexes.map( function(cd_i) { return graph_data1.names_list[cd_i]; });
    data.sort( function(a,b) { return cd_names.indexOf(a.name) - cd_names.indexOf(b.name) });
    draw_plot(data);
});


$("#sort-name").click( function() {
    data.sort( function(a,b) { return a.name.localeCompare(b.name) } );
    draw_plot(data);
});


$("#sort-action").click( function() {
    data.sort( function(a,b) { return (b.cadee_actions.librarymut.length*5 + b.cadee_actions.librarymut_custom.length) - (a.cadee_actions.librarymut.length*5 + a.cadee_actions.librarymut_custom.length) } );
    draw_plot(data);
});


$("#overlay-div").click(function() {
    undock_info_div();
});

// on change for the input textbox - set the cadee action library mut custom value and update the command
$("#CUSTOM_text").on('input propertychange paste', function() {
//    console.log($(this).prop("value"));
    var sel_data = data[x_hovered];
    sel_data.cadee_actions.librarymut_custom = $(this).prop("value");
    update_cadee_cmd();
});

// on change for the librarymut checkboxes - add the cadee action library mut value on the selected mutant
$("#info-actions").find(":checkbox").each(function() {
    $(this).change(function(ev) {
        var sel_data = data[x_hovered];

        if ($(this).prop("id") == "inp_CUSTOM") {
            ct = $("#CUSTOM_text");
            if ($(this).is(":checked")) {
                ct.prop("disabled", false);
                ct.focus();  
                sel_data.cadee_actions.librarymut_custom = ct.prop("value"); 
            } else {
                ct.prop("disabled", true);
                sel_data.cadee_actions.librarymut_custom = "";   // clear the value
            };
        } else {
            j = $.inArray($(this).prop("value"), sel_data.cadee_actions.librarymut);
            if (j > -1) {
                sel_data.cadee_actions.librarymut.splice(j,1);
            } else {
                sel_data.cadee_actions.librarymut.push($(this).prop("value"));
            };
        };
        console.log(sel_data.cadee_actions);
        update_cadee_cmd();
    });
});

// fix for dumbass plotly_click (you have to be really precise to fire the event)
$("#graph").click(function(ev) {
    dock_info_div();
});

$("#graph").mousemove(function(ev) {
    curX = ev.pageX + 20;
    curY = ev.pageY + 20;
    var idiv = $("#info-div");
    if ((curY < 100) || (curY > 450)) {
        idiv.hide();
    } else {
        if (curX > $(window).width()/2) {
            curX -= 540;
        }
        idiv.css({
            "position": "absolute",
            "top": curY + "px",
            "left": curX + "px"
        });    
        if (info_visible) {
            if (!idiv.is(":visible")) {
                idiv.show();
            };
        } else {
            idiv.hide();
        }
    };
});

$(document).keydown(function(e) {
    if (e.keyCode == 16) {
        info_visible = true;
    };    
});

$(document).keyup(function(e){
    if (e.keyCode == 16) {
        info_visible = false;
    }
});

document.getElementById("graph").on('plotly_hover', function(eventData){
    console.log(eventData);
    var tmp = x_hovered;
    x_hovered = Math.round(eventData.xvals[0]);
    if (tmp === x_hovered) { return };
 
    var sel_data = data[x_hovered];
    var cd = document.getElementById("graph").calcdata;  // all three sets
 
    var idiv = $("#info-div");
    $("#info-header").text("System: " + sel_data.name );

    var stats = $("#info-stats");
    var actions = $("#info-actions");

    var html_str = '<p><b>Stats</b></p><table class="pure-table"><thead><tr><th></th><th>Mean</th><th>St.dev.</th><th>Median</th></tr></thead><tbody>';
    ["dG#", "dG0", "dG#_rev"].forEach(function(n,i) {
        var mean = Number(cd[i][x_hovered].mean).toFixed(2);
        var std = Number(cd[i][x_hovered].sd).toFixed(2);
        var med = Number(cd[i][x_hovered].med).toFixed(2);
        html_str = html_str + "<tr><td>" + n + "</td><td>" + mean + "</td><td>" + std + "</td><td>" + med + "</td></tr>";
    });
    html_str = html_str + "</tbody></table>";
    stats.html(html_str);


    actions.find(":checkbox").each(function() {
        i = $.inArray( $(this).prop("value"), sel_data.cadee_actions.librarymut);
        if (i > -1) {
            $(this).prop("checked", true)
        } else {
            $(this).prop("checked", false)
        };

        if ($(this).prop("id") == "inp_CUSTOM") {
            if (sel_data.cadee_actions.librarymut_custom == "") {
                $("#CUSTOM_text").prop("disabled", true);
            } else {
                $("#CUSTOM_text").prop("disabled", false);
                $(this).prop("checked", true)                    
            };
        };
    });
    $("#CUSTOM_text").prop("value", sel_data.cadee_actions.librarymut_custom);    

});

document.getElementById("graph").on('plotly_click', function(eventData){
    console.log(eventData);
    dock_info_div();
});
</script>
"""

if len(sys.argv) != 2:
    print("Usage: \n  " + os.path.basename(__file__) + " cadee.db")
    sys.exit(1)
elif not os.path.lexists(sys.argv[1]):
    print("File %s does not exist!" % sys.argv[1])
    sys.exit(1)


# connect and get values from DB
conn = sqlite3.connect(sys.argv[1])
cursor = conn.cursor()
try:
    cursor.execute("SELECT mutant,barr_forw,exo,barr_back FROM results WHERE feptype='us'")
except sqlite3.DatabaseError as e:
    print("Error when accesing the database: '%s' (%s)" % (sys.argv[1], e))
    sys.exit(1)

results = cursor.fetchall()
conn.close()


# get WT averages
b, e, rb = [], [], []
for res in results:
    mutant, barr, exo, rbarr = res
    if "wt" in mutant.lower():
        b.append(barr)
        e.append(exo)
        rb.append(rbarr)
try:
    AVG_BARR_WT = sum(b)*1.0/len(b)
    AVG_EXO_WT = sum(e)*1.0/len(e)
    AVG_REVBARR_WT = sum(rb)*1.0/len(rb)
except:
    AVG_BARR_WT = 0
    AVG_EXO_WT = 0
    AVG_REVBARR_WT = 0

# get relative energies
data = {}
for res in results:
    mutant, barr, exo, rbarr = res
    mut_name = mutant.split("_")[0].upper()
    if mut_name not in data:
        data[mut_name] = { "name": mut_name, "barrier": [], "exotherm": [], "revbarr": [] }
    barr_rel = round(barr - AVG_BARR_WT, 1)
    exo_rel = round(exo - AVG_EXO_WT, 1)
    rbarr_rel = round(rbarr - AVG_REVBARR_WT, 1)
    data[mut_name]["barrier"].append(barr_rel)
    data[mut_name]["exotherm"].append(exo_rel)
    data[mut_name]["revbarr"].append(rbarr_rel)

dat = """<script> 
var data = {}
</script>
""".format(json.dumps(data.values()))

html += dat + vis + "</body></html>"

open("index.html", 'w').write(html)
print('Success... Wrote index.html... ')
