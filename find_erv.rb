require "./finderv.rb"

if ARGV.length == 1
  puts "Running the following files:"
  Dir[File.join(ARGV[0], "*.xml")].each { |filepath| puts("\t%s" % filepath) }
  
  Finderv.write_putative_erv_batch!(Finderv.putative_ervs_from_batch(Dir[File.join(ARGV[0], "*.xml")]), ARGV[0])
else
  puts "ruby find_erv.rb DIRECTORY_WITH_XML_FILES"
end