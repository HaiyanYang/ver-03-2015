    subroutine initialize_lib_bcd()           
                                              
        integer :: nbcdnodes=0                
        integer :: i=0                        
                                              
     nbcdnodes=324
     allocate(lib_bcdnodes(nbcdnodes))
     lib_bcdnodes=0
     lib_bcdnodes(1)=12
     lib_bcdnodes(2)=13
     lib_bcdnodes(3)=14
     lib_bcdnodes(4)=16
     lib_bcdnodes(5)=19
     lib_bcdnodes(6)=20
     lib_bcdnodes(7)=189
     lib_bcdnodes(8)=190
     lib_bcdnodes(9)=191
     lib_bcdnodes(10)=192
     lib_bcdnodes(11)=193
     lib_bcdnodes(12)=194
     lib_bcdnodes(13)=195
     lib_bcdnodes(14)=196
     lib_bcdnodes(15)=197
     lib_bcdnodes(16)=198
     lib_bcdnodes(17)=199
     lib_bcdnodes(18)=200
     lib_bcdnodes(19)=277
     lib_bcdnodes(20)=278
     lib_bcdnodes(21)=279
     lib_bcdnodes(22)=280
     lib_bcdnodes(23)=281
     lib_bcdnodes(24)=282
     lib_bcdnodes(25)=283
     lib_bcdnodes(26)=284
     lib_bcdnodes(27)=285
     lib_bcdnodes(28)=286
     lib_bcdnodes(29)=287
     lib_bcdnodes(30)=288
     lib_bcdnodes(31)=326
     lib_bcdnodes(32)=327
     lib_bcdnodes(33)=328
     lib_bcdnodes(34)=329
     lib_bcdnodes(35)=330
     lib_bcdnodes(36)=331
     lib_bcdnodes(37)=332
     lib_bcdnodes(38)=333
     lib_bcdnodes(39)=334
     lib_bcdnodes(40)=335
     lib_bcdnodes(41)=336
     lib_bcdnodes(42)=337
     lib_bcdnodes(43)=338
     lib_bcdnodes(44)=339
     lib_bcdnodes(45)=340
     lib_bcdnodes(46)=341
     lib_bcdnodes(47)=342
     lib_bcdnodes(48)=343
     lib_bcdnodes(49)=344
     lib_bcdnodes(50)=345
     lib_bcdnodes(51)=346
     lib_bcdnodes(52)=347
     lib_bcdnodes(53)=348
     lib_bcdnodes(54)=349
     lib_bcdnodes(55)=12
     lib_bcdnodes(56)=13
     lib_bcdnodes(57)=14
     lib_bcdnodes(58)=16
     lib_bcdnodes(59)=19
     lib_bcdnodes(60)=20
     lib_bcdnodes(61)=189
     lib_bcdnodes(62)=190
     lib_bcdnodes(63)=191
     lib_bcdnodes(64)=192
     lib_bcdnodes(65)=193
     lib_bcdnodes(66)=194
     lib_bcdnodes(67)=195
     lib_bcdnodes(68)=196
     lib_bcdnodes(69)=197
     lib_bcdnodes(70)=198
     lib_bcdnodes(71)=199
     lib_bcdnodes(72)=200
     lib_bcdnodes(73)=277
     lib_bcdnodes(74)=278
     lib_bcdnodes(75)=279
     lib_bcdnodes(76)=280
     lib_bcdnodes(77)=281
     lib_bcdnodes(78)=282
     lib_bcdnodes(79)=283
     lib_bcdnodes(80)=284
     lib_bcdnodes(81)=285
     lib_bcdnodes(82)=286
     lib_bcdnodes(83)=287
     lib_bcdnodes(84)=288
     lib_bcdnodes(85)=326
     lib_bcdnodes(86)=327
     lib_bcdnodes(87)=328
     lib_bcdnodes(88)=329
     lib_bcdnodes(89)=330
     lib_bcdnodes(90)=331
     lib_bcdnodes(91)=332
     lib_bcdnodes(92)=333
     lib_bcdnodes(93)=334
     lib_bcdnodes(94)=335
     lib_bcdnodes(95)=336
     lib_bcdnodes(96)=337
     lib_bcdnodes(97)=338
     lib_bcdnodes(98)=339
     lib_bcdnodes(99)=340
     lib_bcdnodes(100)=341
     lib_bcdnodes(101)=342
     lib_bcdnodes(102)=343
     lib_bcdnodes(103)=344
     lib_bcdnodes(104)=345
     lib_bcdnodes(105)=346
     lib_bcdnodes(106)=347
     lib_bcdnodes(107)=348
     lib_bcdnodes(108)=349
     lib_bcdnodes(109)=12
     lib_bcdnodes(110)=13
     lib_bcdnodes(111)=14
     lib_bcdnodes(112)=16
     lib_bcdnodes(113)=19
     lib_bcdnodes(114)=20
     lib_bcdnodes(115)=189
     lib_bcdnodes(116)=190
     lib_bcdnodes(117)=191
     lib_bcdnodes(118)=192
     lib_bcdnodes(119)=193
     lib_bcdnodes(120)=194
     lib_bcdnodes(121)=195
     lib_bcdnodes(122)=196
     lib_bcdnodes(123)=197
     lib_bcdnodes(124)=198
     lib_bcdnodes(125)=199
     lib_bcdnodes(126)=200
     lib_bcdnodes(127)=277
     lib_bcdnodes(128)=278
     lib_bcdnodes(129)=279
     lib_bcdnodes(130)=280
     lib_bcdnodes(131)=281
     lib_bcdnodes(132)=282
     lib_bcdnodes(133)=283
     lib_bcdnodes(134)=284
     lib_bcdnodes(135)=285
     lib_bcdnodes(136)=286
     lib_bcdnodes(137)=287
     lib_bcdnodes(138)=288
     lib_bcdnodes(139)=326
     lib_bcdnodes(140)=327
     lib_bcdnodes(141)=328
     lib_bcdnodes(142)=329
     lib_bcdnodes(143)=330
     lib_bcdnodes(144)=331
     lib_bcdnodes(145)=332
     lib_bcdnodes(146)=333
     lib_bcdnodes(147)=334
     lib_bcdnodes(148)=335
     lib_bcdnodes(149)=336
     lib_bcdnodes(150)=337
     lib_bcdnodes(151)=338
     lib_bcdnodes(152)=339
     lib_bcdnodes(153)=340
     lib_bcdnodes(154)=341
     lib_bcdnodes(155)=342
     lib_bcdnodes(156)=343
     lib_bcdnodes(157)=344
     lib_bcdnodes(158)=345
     lib_bcdnodes(159)=346
     lib_bcdnodes(160)=347
     lib_bcdnodes(161)=348
     lib_bcdnodes(162)=349
     lib_bcdnodes(163)=1
     lib_bcdnodes(164)=2
     lib_bcdnodes(165)=7
     lib_bcdnodes(166)=8
     lib_bcdnodes(167)=23
     lib_bcdnodes(168)=24
     lib_bcdnodes(169)=84
     lib_bcdnodes(170)=85
     lib_bcdnodes(171)=86
     lib_bcdnodes(172)=87
     lib_bcdnodes(173)=88
     lib_bcdnodes(174)=89
     lib_bcdnodes(175)=90
     lib_bcdnodes(176)=91
     lib_bcdnodes(177)=92
     lib_bcdnodes(178)=93
     lib_bcdnodes(179)=94
     lib_bcdnodes(180)=95
     lib_bcdnodes(181)=117
     lib_bcdnodes(182)=118
     lib_bcdnodes(183)=119
     lib_bcdnodes(184)=120
     lib_bcdnodes(185)=121
     lib_bcdnodes(186)=122
     lib_bcdnodes(187)=123
     lib_bcdnodes(188)=124
     lib_bcdnodes(189)=125
     lib_bcdnodes(190)=126
     lib_bcdnodes(191)=127
     lib_bcdnodes(192)=128
     lib_bcdnodes(193)=387
     lib_bcdnodes(194)=388
     lib_bcdnodes(195)=389
     lib_bcdnodes(196)=390
     lib_bcdnodes(197)=391
     lib_bcdnodes(198)=392
     lib_bcdnodes(199)=393
     lib_bcdnodes(200)=394
     lib_bcdnodes(201)=395
     lib_bcdnodes(202)=396
     lib_bcdnodes(203)=397
     lib_bcdnodes(204)=398
     lib_bcdnodes(205)=451
     lib_bcdnodes(206)=452
     lib_bcdnodes(207)=453
     lib_bcdnodes(208)=454
     lib_bcdnodes(209)=455
     lib_bcdnodes(210)=456
     lib_bcdnodes(211)=457
     lib_bcdnodes(212)=458
     lib_bcdnodes(213)=459
     lib_bcdnodes(214)=460
     lib_bcdnodes(215)=461
     lib_bcdnodes(216)=462
     lib_bcdnodes(217)=1
     lib_bcdnodes(218)=2
     lib_bcdnodes(219)=7
     lib_bcdnodes(220)=8
     lib_bcdnodes(221)=23
     lib_bcdnodes(222)=24
     lib_bcdnodes(223)=84
     lib_bcdnodes(224)=85
     lib_bcdnodes(225)=86
     lib_bcdnodes(226)=87
     lib_bcdnodes(227)=88
     lib_bcdnodes(228)=89
     lib_bcdnodes(229)=90
     lib_bcdnodes(230)=91
     lib_bcdnodes(231)=92
     lib_bcdnodes(232)=93
     lib_bcdnodes(233)=94
     lib_bcdnodes(234)=95
     lib_bcdnodes(235)=117
     lib_bcdnodes(236)=118
     lib_bcdnodes(237)=119
     lib_bcdnodes(238)=120
     lib_bcdnodes(239)=121
     lib_bcdnodes(240)=122
     lib_bcdnodes(241)=123
     lib_bcdnodes(242)=124
     lib_bcdnodes(243)=125
     lib_bcdnodes(244)=126
     lib_bcdnodes(245)=127
     lib_bcdnodes(246)=128
     lib_bcdnodes(247)=387
     lib_bcdnodes(248)=388
     lib_bcdnodes(249)=389
     lib_bcdnodes(250)=390
     lib_bcdnodes(251)=391
     lib_bcdnodes(252)=392
     lib_bcdnodes(253)=393
     lib_bcdnodes(254)=394
     lib_bcdnodes(255)=395
     lib_bcdnodes(256)=396
     lib_bcdnodes(257)=397
     lib_bcdnodes(258)=398
     lib_bcdnodes(259)=451
     lib_bcdnodes(260)=452
     lib_bcdnodes(261)=453
     lib_bcdnodes(262)=454
     lib_bcdnodes(263)=455
     lib_bcdnodes(264)=456
     lib_bcdnodes(265)=457
     lib_bcdnodes(266)=458
     lib_bcdnodes(267)=459
     lib_bcdnodes(268)=460
     lib_bcdnodes(269)=461
     lib_bcdnodes(270)=462
     lib_bcdnodes(271)=1
     lib_bcdnodes(272)=2
     lib_bcdnodes(273)=7
     lib_bcdnodes(274)=8
     lib_bcdnodes(275)=23
     lib_bcdnodes(276)=24
     lib_bcdnodes(277)=84
     lib_bcdnodes(278)=85
     lib_bcdnodes(279)=86
     lib_bcdnodes(280)=87
     lib_bcdnodes(281)=88
     lib_bcdnodes(282)=89
     lib_bcdnodes(283)=90
     lib_bcdnodes(284)=91
     lib_bcdnodes(285)=92
     lib_bcdnodes(286)=93
     lib_bcdnodes(287)=94
     lib_bcdnodes(288)=95
     lib_bcdnodes(289)=117
     lib_bcdnodes(290)=118
     lib_bcdnodes(291)=119
     lib_bcdnodes(292)=120
     lib_bcdnodes(293)=121
     lib_bcdnodes(294)=122
     lib_bcdnodes(295)=123
     lib_bcdnodes(296)=124
     lib_bcdnodes(297)=125
     lib_bcdnodes(298)=126
     lib_bcdnodes(299)=127
     lib_bcdnodes(300)=128
     lib_bcdnodes(301)=387
     lib_bcdnodes(302)=388
     lib_bcdnodes(303)=389
     lib_bcdnodes(304)=390
     lib_bcdnodes(305)=391
     lib_bcdnodes(306)=392
     lib_bcdnodes(307)=393
     lib_bcdnodes(308)=394
     lib_bcdnodes(309)=395
     lib_bcdnodes(310)=396
     lib_bcdnodes(311)=397
     lib_bcdnodes(312)=398
     lib_bcdnodes(313)=451
     lib_bcdnodes(314)=452
     lib_bcdnodes(315)=453
     lib_bcdnodes(316)=454
     lib_bcdnodes(317)=455
     lib_bcdnodes(318)=456
     lib_bcdnodes(319)=457
     lib_bcdnodes(320)=458
     lib_bcdnodes(321)=459
     lib_bcdnodes(322)=460
     lib_bcdnodes(323)=461
     lib_bcdnodes(324)=462
    end subroutine initialize_lib_bcd        
