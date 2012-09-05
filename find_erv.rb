require "./finderv.rb"

if ARGV.length == 2 && ARGV[0] =~ /^(erv|ltr)$/i
  puts "Scanning the following files for #{ARGV[0].downcase}:"
  (files = Dir[File.join(ARGV[1], "*.xml")]).each { |filepath| puts("\t%s" % filepath) }
  
  if ARGV[0] =~ /erv/i
    Finderv.write_putative_erv_batch!(Finderv.putative_ervs_from_batch(files), ARGV[1])
  elsif ARGV[0] =~ /ltr/i
    Finderv.write_solo_ltr_batch!(Finderv.solo_ltrs_from_batch(files), ARGV[1])
  end
else
  puts "ruby find_erv.rb (erv OR ltr) DIRECTORY_WITH_XML_FILES"
end
