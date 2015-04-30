# nexustools

Scripts for importing, concatenating, and merging morphological datasets

## TODO

quality control for input datasets

- [x] include character states names
- [x] include character state labels
- [x] make concat script bring in state labels
- [x] identify possible duplicate characters with some symbol to pull out later
- [ ] include character sets for trait types (what are these trait types? maybe take from Livezey and Zusi 2006?)

wishlist

- [ ] Add ability to work with NEXML files
- [ ] Implement concatenate function in R
- [ ] Complete find duplicates function
- [ ] Figure out how to export NEXUS file in R
- [ ] Make a way to scrape text from PDF files?
- [ ] Add ability to subset based on character type
- [ ] Add a way to merge taxa (code overlapping, different scorings as polymorphisms?)
- [ ] Add way to merge duplicate characters (keep most recent? interleave data?)
