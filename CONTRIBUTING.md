## Contributing teaching materials

### Clone Github repository
```bash
# Clone excelerate-scRNAseq
git clone https://github.com/NBISweden/excelerate-scRNAseq.git
cd excelerate-scRNAseq

# Checkout feature branch to work on
git checkout -b session-example
```

### Prepare a self-contained folder with teaching materials for a given session, e.g. session-example
It is probably easiest and most common to use Markdown for the main document, to be able to display content on Github and render it to the website. However, feel free to use other solutions if you have any other strong preferences.

**e.g. folder structure with .md and .Rmd**

```bash
e.g. session-example
.
+-- session-example_files
|   +-- figure-gmf
|       +-- pressure.png
+-- session-example.md
+-- session-example.Rmd
```

*Btw, you get this structure by default if you use R-Studio, New file -> R markdown -> From Template -> GitHub Document (Markdown)*

### Add link(s) to your materials in `schedule.md`
 _see example under Example links_

### Code and commit  
``` bash
 # Code & commit changes while working on the materials
 git add session-feature.md
 git commit -m "commit message"

 # Push to feature when ready
 git push
 ```

_Note Git commit good practices_

**Git commits good practices**
- Commit messages should contain relevant information regarding the feature(s) you add, what type of analyses they can be used for, *etc.*.
- The subject line should be written in an imperative, e.g. *Fix typos* and be 50 characters or less
- The body, if any, should be wrapped at 72 characters.
- The subject and body should be separated by a blank line, start with a capital letter, and the subject line should not end with a period.
- More about [good commit messages][git-commits]

## Go to GitHub to make a pull request to master branch when ready

#### Course website
https://nbisweden.github.io/excelerate-scRNAseq/

#### Questions, feedback etc. 
Create an issue or contact Olga Dethlefsen <olga.dethlefsen@nbis.se>


## [Back to main](README.md)


[git-commits]: https://chris.beams.io/posts/git-commit/
