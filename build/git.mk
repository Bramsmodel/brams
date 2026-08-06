# git info
GIT_BRANCH := $(shell git rev-parse --abbrev-ref HEAD)
GIT_HASH := $(shell git rev-parse --short HEAD)
GIT_DATE := $(shell git log -1 --format=%cd --date=format-local:"%Y-%m-%d @ %H:%M:%S%z")
$(shell echo 'module ModGit' > ${INIT}/ModGit.f90.tmp)
$(shell echo 'character(len=*), parameter :: GitBranch="${GIT_BRANCH}"'>> ${INIT}/ModGit.f90.tmp)
$(shell echo 'character(len=*), parameter :: GitHash="${GIT_HASH}"'>> ${INIT}/ModGit.f90.tmp)
$(shell echo 'character(len=*), parameter :: GitDate="${GIT_DATE}"'>> ${INIT}/ModGit.f90.tmp)
$(shell echo 'end module ModGit' >> ${INIT}/ModGit.f90.tmp)

$(shell \
   echo 'module ModGit' > ${INIT}/ModGit.f90.tmp; \
   echo '  character(len=*), parameter :: GitBranch="${GIT_BRANCH}"' >> ${INIT}/ModGit.f90.tmp; \
   echo '  character(len=*), parameter :: GitHash="${GIT_HASH}"'     >> ${INIT}/ModGit.f90.tmp; \
   echo '  character(len=*), parameter :: GitDate="${GIT_DATE}"'     >> ${INIT}/ModGit.f90.tmp; \
   echo 'end module ModGit' >> ${INIT}/ModGit.f90.tmp; \
   if ! cmp -s ${INIT}/ModGit.f90.tmp ${INIT}/ModGit.f90 2>/dev/null; then \
      mv ${INIT}/ModGit.f90.tmp ${INIT}/ModGit.f90; \
   else \
      rm -f ${INIT}/ModGit.f90.tmp; \
   fi \
)
