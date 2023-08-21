
context('smallestDifferentStrings')

test_that('Common start removed', {
    expect_equal(smallestDifferentStrings(c('START_first', 'START_second', 'START_third'), verbose=F),
                 c('first', 'second', 'third'))
})

test_that('Common end removed', {
    expect_equal(smallestDifferentStrings(c('first_END', 'second_END', 'third_END'), verbose=F),
                 c('first', 'second', 'third'))
})

test_that('Common start and end removed', {
    expect_equal(smallestDifferentStrings(c('START_first_END', 'START_second_END', 'START_third_END'), verbose=F),
                 c('first', 'second', 'third'))
})

