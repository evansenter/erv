#!/usr/bin/env ruby
#
# Sample Ruby script using the Selenium client API
#
require "rubygems"
gem "selenium-client", ">=1.2.16"
require "selenium/client"

begin
  @browser = Selenium::Client::Driver.new \
      :host => "localhost",
      :port => 4444,
      :browser => "*googlechrome",
      :url => "http://www.ncbi.nlm.nih.gov",
      :timeout_in_second => 60

  @browser.start_new_browser_session
  @browser.open "/nuccore/NW_879563.1?report=genbank"
  @browser.click "name=EntrezSystem2.PEntrez.Nuccore.Sequence_ResultsPanel.Sequence_SingleItemSupl.Sequence_ViewerGenbankSidePanel.Sequence_ViewerChangeRegion.Shutter", wait_for: :ajax
  @browser.type "name=EntrezSystem2.PEntrez.Nuccore.Sequence_ResultsPanel.Sequence_SingleItemSupl.Sequence_ViewerGenbankSidePanel.Sequence_ViewerChangeRegion.From", 100
  @browser.type "name=EntrezSystem2.PEntrez.Nuccore.Sequence_ResultsPanel.Sequence_SingleItemSupl.Sequence_ViewerGenbankSidePanel.Sequence_ViewerChangeRegion.To", 1000
  
  sleep 20
  
  # @browser.type "q", "Selenium seleniumhq.org"
  # @browser.click "btnG", :wait_for => :page
  # puts @browser.text?("seleniumhq.org")
ensure
  @browser.close_current_browser_session
end
