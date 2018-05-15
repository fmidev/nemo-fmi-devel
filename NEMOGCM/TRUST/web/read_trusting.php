
<?php

$num=1;

function read_trusting ($file) { 

  global $num;
  $row = 1;

  if (($handle = fopen($file, "r")) !== FALSE) {

    while (($data = fgetcsv($handle, 1000, ";")) !== FALSE) {

      switch ($row) {
        case 1: 
	  $num1 = count($data);
	  print("<thead>\n\n"); print("<tr class=\"title\">\n");
	  break;
        default:	
	  if   ( $row % 2 == 1 ) { print("<tr class=\"odd\">\n") ; }
	  else                   { print("<tr class=\"even\">\n"); }
      }

      $num = count($data);
      #echo "<p> $num champs ligne $row: <br /></p>\n";

      for ($c=0; $c<$num1; $c++) {
        #echo $data[$c] . "<br />\n";
	if ($c > $num ) { $var=""; } else { $var=trim($data[$c]); }

	switch ($var) { 
	  case "FAILED"   : print("<td bgcolor=\"#FF6666\">".$var."</td>\n"); break;
	  case "Different": print("<td bgcolor=\"#FFFF99\">".$var."</td>\n"); break;
	  case "OK"       : print("<td bgcolor=\"#92CA89\">".$var."</td>\n"); break;
	  default         : print("<td>".$var."</td>\n")                           ;
	}

      }

      print("</tr>\n\n");
      if ( $row == 1 ) { print("</thead>\n\n");	print("<tbody>\n\n"); }
      $row++;
    }

    fclose($handle);
  }

  print("</tbody>\n");
}

?>
