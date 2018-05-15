<html>

  <head>

    <title>Trusting Dashboard</title>
    <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
    <script type="text/javascript" src="html/sortabletable.js"></script>
    <link type="text/css" rel="StyleSheet" href="html/sortabletable.css" />
    <link rel="stylesheet" href="html/tab-view.css" type="text/css" media="screen">
    <script type="text/javascript" src="html/tab-view.js"></script>
    <style type="text/css">
      body    { font-family:Helvetica, Arial, Sans-Serif;                      	    	      }
      a1      { color:#000000; text-decoration:none;                           	    	      }
      .title  { background:#AABBBB; font-weight:bold;                          	    	      }
      .even   { background:#eee;                                               	    	      }
      .odd    { background:#ddd;                                               	    	      }
      td      { font-family:Arial    ; font-size:10pt; padding:1px 5px 1px 5px;	    	      }
      .td1    { font-family:Monospace; font-size:10pt; padding:1px 5px 1px 5px;	    	      }
      .table1 { background:#E6EEC9; border-collapse:collapse;                  	    	      }
      .table2 { background:#E6EEC9; border-collapse:collapse;                  	    	      }
      .myh3   { font-weight:bold; font-style:italic; font-family:Arial,Helvetica,sans-serif
                background-color:#FFD700; border:2px outset #FFFFFF
                padding: 1px 20px 1px 20px; width:120px; -moz-border-radius:6px 6px 6px 6px;  }
    </style>

  </head>

  <body>

    <?php
       ini_set('display_errors', 'Off'); error_reporting(E_ALL);
       include 'read_trusting.php';
       $DODS_TGCC="http://dods.extra.cea.fr/store/martin/trusting"; $DODS_IDRIS="http://dodsp.idris.fr/romr005/trusting";
       $ctx=stream_context_create(array('http'=>array('timeout'=>1)));
       $test=file_get_contents("${DODS_TGCC}/ORCA2_LIM_PISCES/nemo_v3_6_STABLE/trusting_info.html", 0, $ctx);
       if ($test == FALSE ) { $LOCATION="./src/martin"; } else { $LOCATION="${DODS_TGCC}"; }
    ?>

    <!--$file_headers=get_headers("${DODS_TGCC}/ORCA2_LIM_PISCES/nemo_v3_6_STABLE/trusting_info.html");
       if ($file_headers[0] == '' ) { $LOCATION="./src/martin"; } else { $LOCATION="${DODS_TGCC}"; }-->

    <!-----Begin NEMO trust----->
    <div id="dhtmlgoodies_tabViewNEMO">

      <!-----Begin v36 trust----->
      <?php $BRANCH="nemo_v3_6_STABLE" ?>
      <div class="dhtmlgoodies_aTab">
	<div id="dhtmlgoodies_tabViewConf">

	  <!-----Begin O2LP trust----->
	  <div class="dhtmlgoodies_aTab">
	    <?php $CONF="ORCA2_LIM_PISCES"; echo file_get_contents("$LOCATION/$CONF/$BRANCH/trusting_info.html"); ?>
	    <div id="dhtmlgoodies_tabViewHPCC">

	      <!-----Begin Curie trust----->
	      <div class="dhtmlgoodies_aTab">
		<table class="sort-table" id="table-1" cellspacing="1">
		  <?php read_trusting("$LOCATION/$CONF/$BRANCH/trusting_martin_X64_CURIE_cron.txt") ?>
		</table>
		<script type="text/javascript">
		  var st = new SortableTable( document.getElementById("table-1"),
	                   [<?php print substr(str_repeat("\"String\",", $num),0,-1); ?>]);
		  st.sort(0, true);
		</script>
	      </div>
	      <!-----End Curie trust----->

	      <!-----Begin Ada trust----->
	      <div class="dhtmlgoodies_aTab">
		<table class="sort-table" id="table-2" cellspacing="1">
		  <?php read_trusting("${DODS_IDRIS}/$CONF/$BRANCH/trusting_romr005_X64_ADA_at.txt") ?>
		</table>
		<script type="text/javascript">
		  var st = new SortableTable( document.getElementById("table-2"),
		           [<?php print substr(str_repeat("\"String\",", $num),0,-1); ?>]);
		  st.sort(0, true);
		</script>
	      </div>
	      <!-----End Ada trust----->

	    </div>
	  </div>
	  <!-----End O2LP trust----->

	  <!-----Begin O1L3P trust----->
	  <div class="dhtmlgoodies_aTab">
	    <?php $CONF="ORCA1_LIM3_PISCES"; echo file_get_contents("$LOCATION/$CONF/$BRANCH/trusting_info.html") ?>
	    <div id="dhtmlgoodies_tabViewHPCC2">

	      <!-----Begin Curie trust----->
	      <div class="dhtmlgoodies_aTab">
		<table class="sort-table" id="table-3" cellspacing="1">
		  <?php read_trusting("$LOCATION/$CONF/$BRANCH/trusting_martin_X64_CURIE_cron.txt") ?>
		</table>
		<script type="text/javascript">
		  var st = new SortableTable( document.getElementById("table-3"),
	                   [<?php print substr(str_repeat("\"String\",", $num),0,-1); ?>]);
		  st.sort(0, true);
		</script>
	      </div>
	      <!-----End Curie trust----->

	      <!-----Begin Ada trust----->
	      <div class="dhtmlgoodies_aTab">
		<table class="sort-table" id="table-4" cellspacing="1">
		  <?php read_trusting("${DODS_IDRIS}/$CONF/$BRANCH/trusting_romr005_X64_ADA_at.txt") ?>
		</table>
		<script type="text/javascript">
		  var st = new SortableTable( document.getElementById("table-4"),
		           [<?php print substr(str_repeat("\"String\",", $num),0,-1); ?>]);
		  st.sort(0, true);
		</script>
	      </div>
	      <!-----End Ada trust----->

	    </div>
	  </div>
	  <!-----End O1L3P trust----->

	  <!-----Begin AMM12 trust----->
	  <div class="dhtmlgoodies_aTab">
	    <?php $CONF="AMM12"; echo file_get_contents("$LOCATION/$CONF/$BRANCH/trusting_info.html") ?>
	    <div id="dhtmlgoodies_tabViewHPCC3">

	      <!-----Begin Curie trust----->
	      <div class="dhtmlgoodies_aTab">
		<table class="sort-table" id="table-5" cellspacing="1">
		  <?php read_trusting("$LOCATION/$CONF/$BRANCH/trusting_martin_X64_CURIE_cron.txt") ?>
		</table>
		<script type="text/javascript">
		  var st = new SortableTable( document.getElementById("table-5"),
	                   [<?php print substr(str_repeat("\"String\",", $num),0,-1); ?>]);
		  st.sort(0, true);
		</script>
	      </div>
	      <!-----End Curie trust----->

	      <!-----Begin Ada trust----->
	      <div class="dhtmlgoodies_aTab">
		<table class="sort-table" id="table-6" cellspacing="1">
		  <?php read_trusting("${DODS_IDRIS}/$CONF/$BRANCH/trusting_romr005_X64_ADA_at.txt") ?>
		</table>
		<script type="text/javascript">
		  var st = new SortableTable( document.getElementById("table-6"),
		           [<?php print substr(str_repeat("\"String\",", $num),0,-1); ?>]);
		  st.sort(0, true);
		</script>
	      </div>
	      <!-----End Ada trust----->

	    </div>
	  </div>
	  <!-----End AMM12 trust----->

	</div>
      </div>
      <!-----End v36 trust----->


      <!-----Begin trunk trust----->
      <?php $BRANCH="trunk" ?>
      <div class="dhtmlgoodies_aTab">
	<div id="dhtmlgoodies_tabViewConf2">

	  <!-----Begin O2LP trust----->
	  <div class="dhtmlgoodies_aTab">
	    <?php $CONF="ORCA2_LIM_PISCES"; echo file_get_contents("$LOCATION/$CONF/$BRANCH/trusting_info.html") ?>
	    <div id="dhtmlgoodies_tabViewHPCC4">

	      <!-----Begin Curie trust----->
	      <div class="dhtmlgoodies_aTab">
		<table class="sort-table" id="table-7" cellspacing="1">
		  <?php read_trusting("$LOCATION/$CONF/$BRANCH/trusting_martin_X64_CURIE_cron.txt") ?>
		</table>
		<script type="text/javascript">
		  var st = new SortableTable( document.getElementById("table-7"),
		           [<?php print substr(str_repeat("\"String\",", $num),0,-1); ?>]);
		  st.sort(0, true);
		</script>
	      </div>
	      <!-----End   Curie trust----->

	      <!-----Begin Ada trust----->
	      <div class="dhtmlgoodies_aTab">
		<table class="sort-table" id="table-8" cellspacing="1">
		  <?php read_trusting("${DODS_IDRIS}/$CONF/$BRANCH/trusting_romr005_X64_ADA_at.txt") ?>
		</table>
		<script type="text/javascript">
		  var st = new SortableTable( document.getElementById("table-8"),
		           [<?php print substr(str_repeat("\"String\",", $num),0,-1); ?>]);
		  st.sort(0, true);
		</script>
	      </div>
	      <!-----End Ada trust----->

	    </div>
	  </div>
	  <!-----End O2LP trust----->

	  <!-----Begin O1L3P trust----->
	  <div class="dhtmlgoodies_aTab">
	    <?php $CONF="ORCA1_LIM3_PISCES"; echo file_get_contents("$LOCATION/$CONF/$BRANCH/trusting_info.html") ?>
	    <div id="dhtmlgoodies_tabViewHPCC5">

	      <!-----Begin Curie trust----->
	      <div class="dhtmlgoodies_aTab">
		<table class="sort-table" id="table-9" cellspacing="1">
		  <?php read_trusting("$LOCATION/$CONF/$BRANCH/trusting_martin_X64_CURIE_cron.txt") ?>
		</table>
		<script type="text/javascript">
		  var st = new SortableTable( document.getElementById("table-9"),
		           [<?php print substr(str_repeat("\"String\",", $num),0,-1); ?>]);
		  st.sort(0, true);
		</script>
	      </div>
	      <!-----End   Curie trust----->

	      <!-----Begin Ada trust----->
	      <div class="dhtmlgoodies_aTab">
		<table class="sort-table" id="table-10" cellspacing="1">
		  <?php read_trusting("${DODS_IDRIS}/$CONF/$BRANCH/trusting_romr005_X64_ADA_at.txt") ?>
		</table>
		<script type="text/javascript">
		  var st = new SortableTable( document.getElementById("table-10"),
		           [<?php print substr(str_repeat("\"String\",", $num),0,-1); ?>]);
		  st.sort(0, true);
		</script>
	      </div>
	      <!-----End Ada trust----->

	    </div>
	  </div>
	  <!-----End O1L3P trust----->

	  <!-----Begin AMM12 trust----->
	  <div class="dhtmlgoodies_aTab">
	    <?php $CONF="AMM12"; echo file_get_contents("$LOCATION/$CONF/$BRANCH/trusting_info.html") ?>
	    <div id="dhtmlgoodies_tabViewHPCC6">

	      <!-----Begin Curie trust----->
	      <div class="dhtmlgoodies_aTab">
		<table class="sort-table" id="table-11" cellspacing="1">
		  <?php read_trusting("$LOCATION/$CONF/$BRANCH/trusting_martin_X64_CURIE_cron.txt") ?>
		</table>
		<script type="text/javascript">
		  var st = new SortableTable( document.getElementById("table-11"),
		           [<?php print substr(str_repeat("\"String\",", $num),0,-1); ?>]);
		  st.sort(0, true);
		</script>
	      </div>
	      <!-----End   Curie trust----->

	      <!-----Begin Ada trust----->
	      <div class="dhtmlgoodies_aTab">
		<table class="sort-table" id="table-12" cellspacing="1">
		  <?php read_trusting("${DODS_IDRIS}/$CONF/$BRANCH/trusting_romr005_X64_ADA_at.txt") ?>
		</table>
		<script type="text/javascript">
		  var st = new SortableTable( document.getElementById("table-12"),
		           [<?php print substr(str_repeat("\"String\",", $num),0,-1); ?>]);
		  st.sort(0, true);
		</script>
	      </div>
	      <!-----End Ada trust----->

	    </div>
	  </div>
	  <!-----End O1L3P trust----->

	</div>
      </div>
      <!-----End trunk trust----->

    </div>
    <!-----End NEMO trust----->

    <!-------------------------------------------------------------------------------------------->
    <script type="text/javascript">
      initTabs('dhtmlgoodies_tabViewNEMO'        ,
               Array('<B>&nbsp;3.6&nbsp;  </B>',
	             '<B>&nbsp;trunk&nbsp;</B>' 
                                                ),
               0,1820,15100                      ,
               Array(false)                       )
      initTabs('dhtmlgoodies_tabViewConf'                                        ,
               Array('<B style="color: green">&nbsp;ORCA2_LIM_PISCES&nbsp;  </B>',
	       	     '<B style="color: green">&nbsp;ORCA1_LIM3_PISCES&nbsp; </B>',
	       	     '<B style="color: green">&nbsp;AMM12&nbsp;             </B>',
	       	     '<B style="color: red  ">&nbsp;C1D_PAPA&nbsp;          </B>',
	       	     '<B style="color: red  ">&nbsp;GYRE&nbsp;              </B>',
	       	     '<B style="color: red  ">&nbsp;GYRE_BFM&nbsp;          </B>',
	       	     '<B style="color: red  ">&nbsp;GYRE_PISCES&nbsp;       </B>',
	       	     '<B style="color: red  ">&nbsp;GYRE_XIOS&nbsp;         </B>',
	       	     '<B style="color: red  ">&nbsp;ORCA2_LIM&nbsp;         </B>',
	       	     '<B style="color: red  ">&nbsp;ORCA2_LIM3&nbsp;        </B>',
	       	     '<B style="color: red  ">&nbsp;ORCA2_LIM_CFC_C14b&nbsp;</B>',
	       	     '<B style="color: red  ">&nbsp;ORCA2_OFF_PISCES&nbsp;  </B>',
	       	     '<B style="color: red  ">&nbsp;ORCA2_SAS_LIM&nbsp;     </B>' 
	                                                                          ),
               0,1810,15050                                                        ,
               Array(false)                                                         )
      initTabs('dhtmlgoodies_tabViewConf2'                                       ,
               Array('<B style="color: green">&nbsp;ORCA2_LIM_PISCES&nbsp;  </B>',
	       	     '<B style="color: green">&nbsp;ORCA1_LIM3_PISCES&nbsp; </B>',
	       	     '<B style="color: green">&nbsp;AMM12&nbsp;             </B>',
	       	     '<B style="color: red  ">&nbsp;C1D_PAPA&nbsp;          </B>',
	       	     '<B style="color: red  ">&nbsp;GYRE&nbsp;              </B>',
	       	     '<B style="color: red  ">&nbsp;GYRE_BFM&nbsp;          </B>',
	       	     '<B style="color: red  ">&nbsp;GYRE_PISCES&nbsp;       </B>',
	       	     '<B style="color: red  ">&nbsp;GYRE_XIOS&nbsp;         </B>',
	       	     '<B style="color: red  ">&nbsp;ORCA2_LIM&nbsp;         </B>',
	       	     '<B style="color: red  ">&nbsp;ORCA2_LIM3&nbsp;        </B>',
	       	     '<B style="color: red  ">&nbsp;ORCA2_LIM_CFC_C14b&nbsp;</B>',
	       	     '<B style="color: red  ">&nbsp;ORCA2_OFF_PISCES&nbsp;  </B>',
	       	     '<B style="color: red  ">&nbsp;ORCA2_SAS_LIM&nbsp;     </B>' 
	                                                                          ),
               0,1810,15050                                                        ,
               Array(false)                                                         )
      initTabs('dhtmlgoodies_tabViewHPCC'                         ,
               Array('<B>&nbsp;TGCC&nbsp;/&nbsp;curie&nbsp;</B>',
	             '<B>&nbsp;IDRIS&nbsp;/&nbsp;ada&nbsp; </B>' ),
               0,1800,15000                                       ,
               Array(false)                                        )
      initTabs('dhtmlgoodies_tabViewHPCC2'                        ,
               Array('<B>&nbsp;TGCC&nbsp;/&nbsp;curie&nbsp;</B>',
	             '<B>&nbsp;IDRIS&nbsp;/&nbsp;ada&nbsp; </B>' ),
               0,1800,15000                                       ,
               Array(false)                                        )
      initTabs('dhtmlgoodies_tabViewHPCC3'                        ,
               Array('<B>&nbsp;TGCC&nbsp;/&nbsp;curie&nbsp;</B>',
	             '<B>&nbsp;IDRIS&nbsp;/&nbsp;ada&nbsp; </B>' ),
               0,1800,15000                                       ,
               Array(false)                                        )
      initTabs('dhtmlgoodies_tabViewHPCC4'                        ,
               Array('<B>&nbsp;TGCC&nbsp;/&nbsp;curie&nbsp;</B>',
	             '<B>&nbsp;IDRIS&nbsp;/&nbsp;ada&nbsp; </B>' ),
               0,1800,15000                                       ,
               Array(false)                                        )
      initTabs('dhtmlgoodies_tabViewHPCC5'                        ,
               Array('<B>&nbsp;TGCC&nbsp;/&nbsp;curie&nbsp;</B>',
	             '<B>&nbsp;IDRIS&nbsp;/&nbsp;ada&nbsp; </B>' ),
               0,1800,15000                                       ,
               Array(false)                                        )
      initTabs('dhtmlgoodies_tabViewHPCC6'                        ,
               Array('<B>&nbsp;TGCC&nbsp;/&nbsp;curie&nbsp;</B>',
	             '<B>&nbsp;IDRIS&nbsp;/&nbsp;ada&nbsp; </B>' ),
               0,1800,15000                                       ,
               Array(false)                                        )
    </script>
    <!-------------------------------------------------------------------------------------------->

  </body>

</html>
