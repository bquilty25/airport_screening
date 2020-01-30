# generate the current date and sink it into an MD file

sprintf("Last built at %s by B. Quilty, S. Clifford, S. Flasche, R. Eggo and other members of CMMID at LSHTM", 
        format(Sys.time(),
               "%d %b %Y at %H:%M:%S")) %>%
  write(., file = "date_stamp.md")


write('<br/><br/>Download the preprint of our analysis <a href="https://github.com/cmmid/cmmid.github.io/raw/master/ncov/airport_screening_report/airport_screening_preprint_2020_01_28.pdf">here</a>.<br>Download the code for this app <a href="https://github.com/bquilty25/airport_screening">on GitHub</a>
<br/><br/>', file = "date_stamp.md", append=T)
