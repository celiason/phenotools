# nexustools

Scripts for importing, concatenating, and merging morphological datasets

## TODO

quality control for input datasets

- [x] include character states names
- [x] include character state labels
- [x] make concat script bring in state labels
- [x] identify possible duplicate characters with some symbol to pull out later

wishlist

- [ ] Make pl script recognize and output character sets from input --> concatenated dataset
- [ ] Add ability to work with NEXML files
- [ ] Implement concatenate function in R
- [ ] Complete find duplicates function
- [ ] Figure out how to export NEXUS file in R
- [ ] Make a way to scrape text from PDF files?
- [ ] Add ability to subset based on character type
- [ ] Add a way to merge taxa (code overlapping, different scorings as polymorphisms?)
- [ ] Add way to merge duplicate characters (keep most recent? interleave data?)
- [ ] Add ability to sort/subset characters based on charsets or grep